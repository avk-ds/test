from picamera2 import Picamera2
from ultralytics import YOLO

# Initialize the Picamera2
picam2 = Picamera2()
picam2.preview_configuration.main.size = (1280, 720)
picam2.preview_configuration.main.format = "RGB888"
picam2.preview_configuration.align()
picam2.configure("preview")
picam2.start()

# Load the YOLO model
model = YOLO("yolo11n.pt")

# Get the center of the frame
frame_width, frame_height = 1280, 720
frame_center = (frame_width // 2, frame_height // 2)

while True:
    # Capture frame
    frame = picam2.capture_array()

    # Run YOLO inference
    results = model(frame)

    # Get object positions
    for result in results:
        for box in result.boxes:
            x1, y1, x2, y2 = map(int, box.xyxy[0])  # Bounding box coordinates
            obj_center = ((x1 + x2) // 2, (y1 + y2) // 2)  # Object center (cx, cy)

            # Print object position
            print(f"Object Center: {obj_center}, Frame Center: {frame_center}")

    # Break loop with 'q'
    if input("Press Enter to continue or 'q' to quit: ") == "q":
        break
