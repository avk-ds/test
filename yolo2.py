import cv2
from picamera2 import Picamera2
from ultralytics import YOLO

# Initialize Picamera2
picam2 = Picamera2()
picam2.preview_configuration.main.size = (1280, 720)
picam2.preview_configuration.main.format = "RGB888"
picam2.preview_configuration.align()
picam2.configure("preview")

# Start the preview
picam2.start()

# Load the YOLO model
model = YOLO("yolo11n.pt")

# Start the camera preview
picam2.start_preview()

while True:
    # Capture frame-by-frame
    frame = picam2.capture_array()

    # Run YOLO inference
    results = model(frame)

    # Get the annotated frame
    annotated_frame = results[0].plot()

    # Update the preview overlay
    picam2.set_overlay(annotated_frame)

