import torch
from . import memory
from . import network
from itertools import count
import random
import numpy
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import math
from datetime import datetime
import rdkit.Chem.Draw
from . import logs
from torch.utils.tensorboard import SummaryWriter
import gc

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class Agent():
    def __init__(self, train_env,
                 num_episodes=5000,  # total number of training episodes
                 replay_memory_length=10000,  # max num transitions in buffer
                 backprop_period=1,  # number of episodes between backprops
                 batch_size=128,  # transitions per SGD batch
                 discount_factor=0.999,  # contribution of future rewards to Q
                 learning_rate=0.0001,  # fraction of gradient to add to params
                 tau=0.001,  # fraction of gradient update to apply to target
                 logging=True):  # save logging plots
        self.num_episodes = num_episodes
        self.train_env = train_env
        # memories
        self.memory = memory.ReplayMemory(replay_memory_length)
        # networks
        self.batch_size = batch_size
        self.policy_net = network.FingerprintNetwork().to(device)
        self.target_network = network.FingerprintNetwork().to(device)
        self.target_network.load_state_dict(self.policy_net.state_dict())
        self.target_network.eval()
        # network optimization
        self.gamma = discount_factor
        self.optimizer_lr = learning_rate
        self.optimizer = torch.optim.Adam(self.policy_net.parameters(),
                                          lr=self.optimizer_lr)
        self.tau = tau  # soft update interpolation parameter
        self.backprop_period = backprop_period
        # logging
        self.writer = SummaryWriter()
        self.writer.add_graph(self.policy_net, torch.rand(2049, device=device))
        self.logging_level = 2  # options 1: minimal figures, 2: some, 3: all
        self.eval_period = 10
        self.plot_period = 50
        self.logging = logging
        self.save_model = False
        self.save_period = 100
        self.num_Qs = 0
        # epsilon decay parameters
        self.initial_epsilon = 1
        self.final_epsilon = 0.01
        self.decay_endpoint = int(0.25*self.num_episodes)
        self.epsilon_gradient = ((self.final_epsilon - self.initial_epsilon) /
                                 self.decay_endpoint)
        self.epsilon_decay = num_episodes/5

    def print_network(self):
        print(self.policy_net)

    def select_action(self, options, epsilon):
        if random.random() > epsilon:
            with torch.no_grad():
                Q = self.policy_net(options.to(device))
                if self.logging and self.logging_level >= 2:
                    self.num_Qs += 1
                    self.writer.add_scalar('Q* per Non-Stochastic Action',
                                           Q.max(0)[0].detach(),
                                           self.num_Qs)
                return Q.max(0)[1].item()
        else:
            return random.randrange(options.size()[0])

    def optimize_model(self):
        if len(self.memory) < self.batch_size:
            return
        transitions = self.memory.sample(self.batch_size)
        batch = memory.Transition(*zip(*transitions))

        q_tp1 = self.policy_net(torch.stack(batch.next_state).to(device))
        max_q_tp1_prime = torch.zeros(
            self.batch_size, 1, requires_grad=False).to(device)
        for i in range(self.batch_size):
            max_q_tp1_prime[i] = self.target_network(
                batch.next_state_options[i].to(device)
            ).max(0)[0].detach()

        max_q_tp1_prime_masked = max_q_tp1_prime * \
            torch.cat(batch.terminal).to(device)
        y = torch.stack(batch.reward).to(device).add(
            self.gamma*max_q_tp1_prime_masked)

        loss = torch.nn.functional.huber_loss(q_tp1, y)

        self.optimizer.zero_grad(set_to_none=True)
        loss.backward()
        self.optimizer.step()
        self.soft_update()
        for p in self.policy_net.parameters():
            p.grad = None
        del transitions
        del batch
        del q_tp1
        del max_q_tp1_prime
        del max_q_tp1_prime_masked
        del y
        return loss.detach()

    def soft_update(self):
        with torch.no_grad():
            for target_param, policy_param in zip(self.target_network.parameters(),
                                                  self.policy_net.parameters()):
                target_param.data.copy_(self.tau*policy_param.data +
                                        (1-self.tau)*target_param.data)

    def shape_step_reward(self, final_state_reward, steps_remaining):
        rf = final_state_reward
        tmax = self.train_env.config['training']['max_steps']
        ti = steps_remaining
        return (rf - (rf/tmax)*(ti-1))

    def shape_final_reward(self, final_state_reward):
        rf = final_state_reward
        if rf >= 12:
            return rf
        else:
            return 6*rf - 60

    def update_epsilon_linear(self, i):
        return self.epsilon_gradient*i + self.initial_epsilon

    def update_epsilon_exponential(self, i):
        return math.exp(-1.*(i/self.epsilon_decay))

    def evaluate(self, episode):
        options = self.train_env.reset()
        for t in count():
            action = self.select_action(options, epsilon=0)
            options, reward, terminal = self.train_env.step(action)
            if terminal == 0:
                self.writer.add_scalar(
                    'Terminal Evaluation BA v. Training Episode',
                    reward, episode)
                if reward >= 13:
                    self.writer.add_image(
                        'Eval Mol {}'.format(episode),
                        self.train_env.get_state_image_tensor())
                break

    def train(self):
        # clear terminal
        os.system('cls' if os.name == 'nt' else 'clear')
        # iterate through episodes
        for i_episode in tqdm(range(self.num_episodes), desc="Training"):
            # initialize temp memory buffer
            temp_memory = []
            # update epsilon
            if i_episode <= self.decay_endpoint:
                epsilon = self.update_epsilon_linear(i_episode)
            else:
                epsilon = self.final_epsilon
            # initialize reward total for dense rewards
            episode_total_reward = torch.tensor(
                [0], dtype=torch.float32, device=device)
            # initialize environment
            options = self.train_env.reset()
            # save the starting state and options
            if self.logging and self.logging_level == 3:
                self.writer.add_image(
                    "Train Episode {} State {}".format(i_episode+1, 0),
                    self.train_env.get_state_image_tensor())
                self.writer.add_image(
                    "Train Episode {} State {} Options".format(i_episode+1, 0),
                    self.train_env.get_options_image_tensor())
            # episode loop
            for t in count():
                # act and receive feedback from environment
                action = self.select_action(options, epsilon)
                next_state = options[action]
                options, reward, done = self.train_env.step(action)
                # save the current state and options
                if self.logging and self.logging_level == 3:
                    self.writer.add_image(
                        "Train Episode {} State {}".format(i_episode+1, t+1),
                        self.train_env.get_state_image_tensor())
                    self.writer.add_image(
                        "Train Episode {} State {} Options".format(
                            i_episode+1, t+1),
                        self.train_env.get_options_image_tensor())
                # save transition to temp memory buffer
                temp_memory.append(
                    memory.Transition(
                        next_state,  # torch.tensor
                        options,  # torch.tensor((torch.tensor))
                        reward,  # torch.tensor([float])
                        done,  # torch.tensor(boolean)
                    )
                )
                # if the last transition was terminal
                if done == 0:
                    for t in temp_memory:
                        step_reward = self.shape_step_reward(reward,
                                                             t.next_state[-1])
                        self.memory.push(t.next_state,
                                         t.next_state_options,
                                         step_reward,
                                         t.terminal)
                        episode_total_reward += step_reward
                        del t
                    if (i_episode % self.backprop_period == 0):
                        loss = self.optimize_model()
                    if self.save_model and (i_episode % self.save_period == 0):
                        torch.save(self.policy_net.state_dict(),
                                   "policy-parameters.pt")
                    if self.logging:
                        if loss:  # if batch_size > len(memory) loss = None
                            self.writer.add_scalar(
                                "Mean Batch Loss v. Training Episode",
                                loss, i_episode)
                        if self.logging_level >= 2:
                            self.writer.add_scalar(
                                "Epsilon v. Training Episode",
                                epsilon, i_episode)
                            self.writer.add_scalar(
                                "Terminal Training BA v. Training Episode",
                                reward, i_episode)
                            self.writer.add_scalar(
                                "Total Training Episode Reward v. Training Episode",
                                episode_total_reward, i_episode)
                            if reward >= 13:
                                self.writer.add_image(
                                    'Train Mol {}'.format(i_episode),
                                    self.train_env.get_state_image_tensor())
                        if i_episode % self.eval_period == 0:
                            # only evaluate after saving the final env state
                            # as it resets the environment
                            self.evaluate(i_episode)
                    os.system('rm ligand_*')
                    os.system('rm *map*')
                    os.system('rm *gpf*')
                    del episode_total_reward
                    del loss
                    del action
                    del next_state
                    del options
                    del reward
                    del done
                    del temp_memory
                    gc.collect()
                    torch.cuda.empty_cache()
                    break

        print('\nComplete')
