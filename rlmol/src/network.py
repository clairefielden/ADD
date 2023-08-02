import torch


class FingerprintNetwork(torch.nn.Module):
    def __init__(self):
        super(FingerprintNetwork, self).__init__()
        self.linear1 = torch.nn.Linear(2048+1, 1024)
        self.linear2 = torch.nn.Linear(1024, 512)
        self.linear3 = torch.nn.Linear(512, 128)
        self.linear4 = torch.nn.Linear(128, 32)
        self.linear5 = torch.nn.Linear(32, 1)

    def forward(self, x):
        # x: batch of states accessible from the current state, each
        # concatenated with the steps remaining from the current state.
        h_relu = torch.nn.functional.relu(self.linear1(x))
        h_relu = torch.nn.functional.relu(self.linear2(h_relu))
        h_relu = torch.nn.functional.relu(self.linear3(h_relu))
        h_relu = torch.nn.functional.relu(self.linear4(h_relu))
        h_relu = self.linear5(h_relu)
        return h_relu
