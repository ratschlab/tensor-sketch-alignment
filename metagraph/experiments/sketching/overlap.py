import plotly.graph_objects as go
import json
from plotly.subplots import make_subplots

if __name__ == '__main__':
    K_VALS = list(range(11, 90, 10))
    
    SKETCH_FILE = "./runs/2022-05-21T01:09:23/points.json"
    EXACT_FILE = "./runs/2022-05-21T01:02:45/points.json"
    
    sketch_data = json.load(open(SKETCH_FILE))
    exact_data = json.load(open(EXACT_FILE))
    fig = make_subplots(rows=2, cols=1)
    fig.add_trace(go.Scatter(x=sketch_data['avg_time'], y=sketch_data['recall'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Sketch"), row=1, col=1)

    fig.add_trace(go.Scatter(x=exact_data['avg_time'], y=exact_data['recall'], text=K_VALS, textposition="top center", mode="lines+markers+text", name="Harun"), row=2, col=1)

    fig.update_layout(
        title = "Seed recall",
        xaxis_title = "Avg time",
        yaxis_title = "Recall"
        )

    fig.write_image("overlap.png", scale=1, width=1920, height=1080)
