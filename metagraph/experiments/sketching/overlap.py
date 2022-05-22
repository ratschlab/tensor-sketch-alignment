import plotly.graph_objects as go
import json
from plotly.subplots import make_subplots
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-k", type=int, required=True)
    parser.add_argument("--sketch-json", type=str, required=True)
    parser.add_argument("--exact-json", type=str, required=True)
    args = parser.parse_args()

    MAX_K = args.max_k
    SKETCH_FILE = args.sketch_json
    EXACT_FILE = args.exact_json

    K_VALS = list(range(11, MAX_K, 10))
    
    sketch_data = json.load(open(SKETCH_FILE))
    exact_data = json.load(open(EXACT_FILE))

    fig = make_subplots(rows=2, cols=1, subplot_titles=["Recall for sketch seeder @ varying k", "Recall for default seeder @ varying k"])
    fig.add_trace(go.Scatter(x=sketch_data['avg_time'], y=sketch_data['recall'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Sketch"), row=1, col=1)
    fig.add_trace(go.Scatter(x=exact_data['avg_time'], y=exact_data['recall'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Default"), row=2, col=1)
    fig.update_xaxes(title_text="Average Time (s)", row=1, col=1)
    fig.update_xaxes(title_text="Average Time (s)", row=2, col=1)
    fig.update_yaxes(title_text="Recall", row=1, col=1)
    fig.update_yaxes(title_text="Recall", row=2, col=1)
    fig.write_image("recall.png", scale=1, width=1920, height=1080)

    fig = make_subplots(rows=2, cols=1, subplot_titles=["Fwd precision for sketch seeder @ varying k", "Fwd precision for default seeder @ varying k"])
    fig.add_trace(go.Scatter(x=sketch_data['avg_time'], y=sketch_data['fwd_precision'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Sketch"), row=1, col=1)
    fig.add_trace(go.Scatter(x=exact_data['avg_time'], y=exact_data['fwd_precision'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Default"), row=2, col=1)
    fig.update_xaxes(title_text="Average Time (s)", row=1, col=1)
    fig.update_xaxes(title_text="Average Time (s)", row=2, col=1)
    fig.update_yaxes(title_text="Fwd Precision", row=1, col=1)
    fig.update_yaxes(title_text="Fwd Precision", row=2, col=1)
    fig.write_image("fwd_precision.png", scale=1, width=1920, height=1080)

    fig = make_subplots(rows=2, cols=1, subplot_titles=["RC precision for sketch seeder @ varying k", "RC precision for default seeder @ varying k"])
    fig.add_trace(go.Scatter(x=sketch_data['avg_time'], y=sketch_data['rc_precision'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Sketch"), row=1, col=1)
    fig.add_trace(go.Scatter(x=exact_data['avg_time'], y=exact_data['rc_precision'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Default"), row=2, col=1)
    fig.update_xaxes(title_text="Average Time (s)", row=1, col=1)
    fig.update_xaxes(title_text="Average Time (s)", row=2, col=1)
    fig.update_yaxes(title_text="RC Precision", row=1, col=1)
    fig.update_yaxes(title_text="RC Precision", row=2, col=1)
    fig.write_image("rc_precision.png", scale=1, width=1920, height=1080)
