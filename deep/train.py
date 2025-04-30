# # train.py
#
# import torch
# from torch import nn, optim
# from models.sramp import SRAMP
# from utils.data_loader import load_data
# from config import params
#
# for i in range(1, 11):
#     print(f"🔁 Training run {i}")
#     x_train, x_test, emb_train, emb_test, y_train, y_test, train_loader = load_data(i)
#
#     model = SRAMP(mode='genomeonly')
#     optimizer = optim.Adam(model.parameters(), lr=params.learning_rate)
#     criterion = nn.BCELoss()
#
#     for epoch in range(params.epochs):
#         model.train()
#         for batch_x, batch_emb, batch_y in train_loader:
#             pred = model(batch_x, batch_emb)
#             loss = criterion(pred, batch_y)
#             optimizer.zero_grad()
#             loss.backward()
#             optimizer.step()
#
#     torch.save(model.state_dict(), f"sramp_run{i}.pt")
#     print(f"✅ Saved model: sramp_run{i}.pt")
import torch
from torch import nn, optim
from config.params import Config
from models.sramp import SRAMP
from utils.data_loader import load_data, get_loader
from utils.metrics import evaluate
import pandas as pd

summary_results = {
    "round": [],
    "mode": [],
    "accuracy": [],
    "recall": [],
    "mcc": [],
    "f1": [],
    "auc": [],
}

for i in range(1, 11):
    x_train, x_test, emb_train, emb_test, y_train, y_test = load_data(i)

    for mode in ['genomeonly', 'seqonly']:
        print(f"🔁 Training round {i} | mode: {mode}")
        train_loader = get_loader(x_train, emb_train, y_train, Config.batch_size)

        model = SRAMP(mode=mode)
        criterion = nn.BCELoss()
        optimizer = optim.Adam(model.parameters(), lr=Config.lr)

        for epoch in range(Config.epochs):
            model.train()
            for bx, be, by in train_loader:
                pred = model(bx, be)
                loss = criterion(pred, by)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

        # Save model
        model_path = f'model_run{i}_{mode}.pt'
        torch.save(model.state_dict(), model_path)
        print(f"✅ Saved {model_path}")

        # Evaluation
        model.eval()
        with torch.no_grad():
            y_pred = model(x_test, emb_test).cpu().numpy()
            y_true = y_test.cpu().numpy()
            metrics = evaluate(y_true, y_pred)

        # Print and save metrics
        print(f"📊 Eval round {i} | mode: {mode} | {metrics}")
        summary_results["round"].append(i)
        summary_results["mode"].append(mode)
        for k in metrics:
            summary_results[k].append(metrics[k])

# Create and print summary DataFrame
df = pd.DataFrame(summary_results)
summary_df = df.groupby("mode").mean(numeric_only=True)
print("\n🔍 Summary Metrics (averaged over 10 runs):")
print(summary_df.round(4))

# Save summary if needed
df.to_csv("all_metrics_by_round.csv", index=False)
summary_df.to_csv("summary_by_mode.csv")
