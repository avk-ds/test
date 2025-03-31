import cv2
import numpy as np
import matplotlib.pyplot as plt
from picamera2 import Picamera2
from ultralytics import YOLO

# Initialize the Picamera2
picam2 = Picamera2()
picam2.preview_configuration.main.size = (1280, 720)
picam2.preview_configuration.main.format = "RGB888"
picam2.preview_configuration.align()
picam2.configure("preview")
picam2.start()

# Load the YOLO11 model
model = YOLO("yolo11n.pt")

while True:
    # Capture frame-by-frame
    frame = picam2.capture_array()

    # Run YOLO11 inference on the frame
    results = model(frame)

    # Visualize the results on the frame
    annotated_frame = results[0].plot()

    # Convert to RGB for matplotlib
    annotated_frame = cv2.cvtColor(annotated_frame, cv2.COLOR_BGR2RGB)

    # Display with matplotlib (since cv2.imshow doesn't work on Raspberry Pi OS Lite)
    plt.imshow(annotated_frame)
    plt.axis("off")
    plt.pause(0.01)  # Small pause to refresh the image

    # Exit condition (press Ctrl+C to stop)
