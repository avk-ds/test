import pygame
import numpy as np
from picamera2 import Picamera2
from ultralytics import YOLO

# Initialize Pygame
pygame.init()
screen_size = (1280, 720)
screen = pygame.display.set_mode(screen_size)
pygame.display.set_caption("YOLO Object Tracking")

# Initialize Picamera2
picam2 = Picamera2()
picam2.preview_configuration.main.size = screen_size
picam2.preview_configuration.main.format = "RGB888"
picam2.preview_configuration.align()
picam2.configure("preview")
picam2.start()

# Load YOLO model
model = YOLO("yolo11n.pt")

running = True
while running:
    # Capture frame
    frame = picam2.capture_array()

    # Run YOLO inference
    results = model(frame)
    annotated_frame = results[0].plot()

    # Convert the annotated frame to a format suitable for Pygame
    annotated_frame = np.rot90(annotated_frame)  # Rotate if needed
    annotated_frame = np.flip(annotated_frame, axis=0)  # Flip to match display
    annotated_frame = pygame.surfarray.make_surface(annotated_frame)

    # Display the frame in Pygame window
    screen.blit(annotated_frame, (0, 0))
    pygame.display.flip()

    # Check for quit event
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

# Clean up
pygame.quit()
