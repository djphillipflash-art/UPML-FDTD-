import torch
print("torch版本", torch.__version__)
print("mps可用", torch.backends.mps.is_available())