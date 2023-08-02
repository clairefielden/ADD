import torch
import numpy
import matplotlib.pyplot as plt
import os


def log_q_value(q_value_log):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(q_value_log) + 1, 1)
    ax.plot(x_range, q_value_log, 'mo--')
    ax.set(xlabel='Non-stochastic Step', ylabel='Q-Value',
           title='Q-Values when the agent acts non-stochastically')

    fig.savefig(os.path.join(".", "training-figures",
                "q-value-log.png"))
    plt.close(fig)


def plot_durations(episode_durations):
    fig, ax = plt.subplots()
    durations_t = torch.tensor(episode_durations, dtype=torch.float)
    x_range = numpy.arange(1, len(episode_durations) + 1, 1)
    ax.plot(x_range, durations_t.numpy(), 'mo--')
    ax.set(xlabel='Training Episode', ylabel='Steps',
           title='Steps Per Training Episode')
    # ax.set_xticks(x_range)

    if len(durations_t) >= 100:
        means = durations_t.unfold(0, 100, 1).mean(1).view(-1)
        means = torch.cat((torch.zeros(99), means))
        ax.plot(means.numpy())

    fig.savefig(os.path.join(".", "training-figures",
                "training-episode-durations.png"))
    plt.close(fig)


def plot_total_reward(episode_total_rewards):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(episode_total_rewards) + 1, 1)
    ax.plot(x_range, episode_total_rewards, 'mo--')
    ax.set(xlabel='Training Episode', ylabel='Total Reward',
           title='Total Reward Per Training Episode')
    # ax.set_xticks(x_range)

    fig.savefig(os.path.join(".", "training-figures",
                "training-episode-total-rewards.png"))
    plt.close(fig)


def plot_binding_affinities(binding_affinities):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(binding_affinities) + 1, 1)
    ax.plot(x_range, binding_affinities, 'mo--')
    ax.set(xlabel='Training Episode', ylabel='Binding Affinity (kcal/mol)',
           title='Binding Affinity of Final State per Training Episode')
    # ax.set_xticks(x_range)

    fig.savefig(os.path.join(".", "training-figures",
                "binding_affinities.png"))
    plt.close(fig)


def plot_QEDs(qed_list):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(qed_list) + 1, 1)
    ax.plot(x_range, qed_list, 'mo--')
    ax.set(xlabel='Training Episode', ylabel='QED',
           title='QED of Final State per Training Episode')
    # ax.set_xticks(x_range)

    fig.savefig(os.path.join(".", "training-figures",
                "training-episode-final-qed.png"))
    plt.close(fig)


def plot_eval_QEDs(qed_list):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(qed_list) + 1, 1)
    ax.plot(x_range, qed_list, 'mo--')
    ax.set(xlabel='Evaluation Episode', ylabel='QED',
           title='QED of Final State per Evaluation Episode')
    # ax.set_xticks(x_range)

    fig.savefig(os.path.join(".", "training-figures",
                "eval-episode-final-qed.png"))
    plt.close(fig)


def plot_avg_eval_rewards(avg_eval_rewards, eval_period, num_eval_episodes):
    fig, ax = plt.subplots()
    x_range = numpy.arange(eval_period,
                           (eval_period * len(avg_eval_rewards)
                            + eval_period),
                           eval_period)
    ax.plot(x_range, avg_eval_rewards, 'mo--')
    ax.set(xlabel='Training Episode',
           ylabel='Average Evaluation Environment Reward',
           title="Average Evaluation Reward Per " +
           "{} Evaluation Episodes".format(num_eval_episodes))

    fig.savefig(os.path.join(".", "training-figures",
                "average-evaluation-reward.png"))
    plt.close(fig)


def plot_losses(losses, optimizer_lr):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(losses) + 1, 1)
    ax.plot(x_range, losses, 'mo--')
    ax.set(xlabel='Backprop Step', ylabel='Loss',
           title='Loss at Backprop Step ' +
           '(lr Adam = {}) '.format(optimizer_lr))
    fig.savefig(os.path.join(".", "training-figures",
                "loss.png"))
    plt.close(fig)


def plot_epsilon(epsilons, epsilon_decay):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(epsilons) + 1, 1)
    ax.plot(x_range, epsilons, 'mo--')
    ax.set(xlabel='Training Episode', ylabel='Epsilon',
           title='Stochastic Probability per Training Episode (Decay = {})'
           .format(epsilon_decay))
    fig.savefig(os.path.join(".", "training-figures",
                "training-episode-epsilon.png"))
    plt.close(fig)


def plot_latent_distances(distances):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(distances) + 1, 1)
    ax.plot(x_range, distances, 'mo--')
    ax.set(xlabel='Training Episode', ylabel='Latent Space Distance',
           title='Latent Space Distance at end of Training Episode')
    # ax.set_ylim([0, self.train_env.goal_radius])
    fig.savefig(os.path.join(".", "training-figures",
                "training-episode-latent-space-distance.png"))
    plt.close(fig)


def plot_train_similarities(train_similarities):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(train_similarities) + 1, 1)
    ax.plot(x_range, train_similarities, 'mo--')
    ax.set(xlabel='Training Episode', ylabel='Tanimoto Similarity',
           title='Tanimoto Similarity at end of Training Episode')
    # ax.set_ylim([0, self.train_env.goal_radius])
    fig.savefig(os.path.join(".", "training-figures",
                "training-episode-tanimoto-similarity.png"))
    plt.close(fig)


def plot_eval_distances(eval_distances):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(eval_distances) + 1, 1)
    ax.plot(x_range, eval_distances, 'mo--')
    ax.set(xlabel='Eval Episode', ylabel='Latent Space Distance',
           title='Latent Space Distance at end of ' +
           'Evaluation Episode')
    fig.savefig(os.path.join(".", "training-figures",
                "eval-episode-latent-space-distance.png"))
    plt.close(fig)


def plot_eval_similarities(eval_similarities):
    fig, ax = plt.subplots()
    x_range = numpy.arange(1, len(eval_similarities) + 1, 1)
    ax.plot(x_range, eval_similarities, 'mo--')
    ax.set(xlabel='Eval Episode', ylabel='Tanimoto Similarity',
           title='Tanimoto Similarity at end of ' +
           'Evaluation Episode')
    fig.savefig(os.path.join(".", "training-figures",
                "eval-episode-tanimoto-similarity.png"))
    plt.close(fig)


def save_final_mol(train_env, episode):
    train_env.save(os.path.join(".", "training-figures",
                   "training_episode_{}.png".format(episode + 1)))


def save_final_mol_in_range(train_env, episode):
    train_env.save_previous(os.path.join(".", "training-figures",
                            "training_episode_{}.png".format(episode + 1)))


def save_final_eval_mol(eval_env, episode):
    eval_env.save(os.path.join(".", "training-figures",
                               "evaluation_episode_{}.png".format(episode + 1)))


def log_state(env, path, t):
    # save the current state and options
    env.save(os.path.join(path, 'test_step_{}_state.png'.format(t+1)))


def log_state_options(env, path, t):
    env.saveStateOptions(os.path.join(path,
                         'test_step_{}_state_options.png'.format(t+1)))


def log_starting_QED(env, dir_path, i, start_qed):
    path = os.path.join(
        dir_path, "mol_{}_a_start-qed_{:.2f}.png".format(i, start_qed))
    env.save(path)


def log_ending_QED(env, dir_path, i, end_qed):
    path = os.path.join(
        dir_path, "mol_{}_b_end-qed_{:.2f}.png".format(i, end_qed))
    env.save(path)
