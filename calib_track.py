import time
import numpy as np
from picamera2 import Picamera2
from ultralytics import YOLO
from pyfirmata2 import Arduino, SERVO

# Initialize the Arduino
arduino = Arduino('/dev/ttyUSB0')  # Change to your actual port
servo_pin = 9  # Define the servo control pin
laser_pin = 8  # Define the laser control pin

# Set pin modes
arduino.digital[servo_pin].mode = SERVO
arduino.digital[laser_pin].mode = 1  # Set laser pin as OUTPUT

# Initialize the Picamera2
picam2 = Picamera2()
picam2.preview_configuration.main.size = (1280, 720)
picam2.preview_configuration.main.format = "RGB888"
picam2.preview_configuration.align()
picam2.configure("preview")
picam2.start()

# Load the YOLO model
model = YOLO("yolo11n.pt")

# Frame dimensions
frame_width, frame_height = 1280, 720
frame_center = frame_width // 2  # Only considering X-axis for alignment

# Servo initial position
servo_angle = 90  # Start at the middle position
arduino.digital[servo_pin].write(servo_angle)

# Calibration data storage
frame_x_positions = []
servo_angles = []

def wait_for_target():
    """Waits until an object is detected and returns its X-coordinate."""
    print("Waiting for target detection...")
    while True:
        frame = picam2.capture_array()
        results = model(frame)

        for result in results:
            for box in result.boxes:
                x1, _, x2, _ = map(int, box.xyxy[0])  # Get bounding box
                obj_x = (x1 + x2) // 2  # Compute object center (only X-axis)
                print(f"Target detected at X: {obj_x}")
                return obj_x  # Return the detected object's X position

def calibrate():
    """Perform calibration by detecting targets and allowing trial-and-error angle adjustments."""
    global servo_angle
    num_points = 5  # Number of calibration points

    for i in range(num_points):
        obj_x = wait_for_target()  # Wait for object detection
        arduino.digital[laser_pin].write(1)  # Turn laser ON
        time.sleep(1)

        print(f"Laser fired from default servo angle ({servo_angle}°).")
        print("Keep entering angles until the laser is perfectly aligned.")

        while True:
            try:
                new_angle = float(input(f"Enter servo angle (current: {servo_angle}°): "))
                if 0 <= new_angle <= 180:  # Ensure within servo range
                    servo_angle = new_angle
                    arduino.digital[servo_pin].write(servo_angle)  # Move servo
                    print(f"Servo moved to {servo_angle}°")
                else:
                    print("Angle must be between 0° and 180°.")
                
                confirm = input("Is the laser aligned? (y/n): ").strip().lower()
                if confirm == 'y':
                    break  # Exit the trial-and-error loop when satisfied
                
            except ValueError:
                print("Invalid input. Please enter a number.")

        # Store calibration data
        frame_x_positions.append(obj_x)
        servo_angles.append(servo_angle)

        arduino.digital[laser_pin].write(0)  # Turn laser OFF
        time.sleep(1)

    # Fit a linear function: Servo Angle = m * Object X + b
    m, b = np.polyfit(frame_x_positions, servo_angles, 1)
    
    print("\nCalibration Complete!")
    print(f"Calibration equation: Servo Angle = {m:.4f} * Object X + {b:.4f}")

    return m, b

def track_objects(m, b):
    """Tracks objects and moves the servo to align the camera with the detected object."""
    print("\nTracking started...")
    while True:
        frame = picam2.capture_array()
        results = model(frame)

        for result in results:
            for box in result.boxes:
                x1, _, x2, _ = map(int, box.xyxy[0])  # Get bounding box
                obj_x = (x1 + x2) // 2  # Compute object center (only X-axis)
                
                # Compute the required servo angle using calibration data
                new_servo_angle = int(m * obj_x + b)
                new_servo_angle = max(0, min(180, new_servo_angle))  # Ensure within servo limits

                print(f"Object at X: {obj_x}, Moving servo to {new_servo_angle}°")
                arduino.digital[servo_pin].write(new_servo_angle)

        if input("Press Enter to continue or 'q' to quit: ") == "q":
            break

# Run Calibration
m, b = calibrate()

# Start Tracking
track_objects(m, b)
