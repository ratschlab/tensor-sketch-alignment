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

    K_VALS = list(range(10, MAX_K, 10))
    
    sketch_data = json.load(open(SKETCH_FILE))
    exact_data = json.load(open(EXACT_FILE))
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sketch_data['avg_time'], y=sketch_data['recall'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Sketch"))
    fig.add_trace(go.Scatter(x=exact_data['avg_time'], y=exact_data['recall'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Default"))
    fig.update_xaxes(title_text="Average Time (s)", type="log")
    fig.update_yaxes(title_text="Recall")
    fig.write_image("recall.png", scale=1, width=1920, height=1080)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sketch_data['avg_time'], y=sketch_data['precision'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Sketch"))
    fig.add_trace(go.Scatter(x=exact_data['avg_time'], y=exact_data['precision'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Default"))
    fig.update_xaxes(title_text="Average Time (s)", type="log")
    fig.update_yaxes(title_text="Average #kmers/node")
    fig.write_image("precision.png", scale=1, width=1920, height=1080)
