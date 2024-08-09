import serial
import socket
import sys
import threading

# Serial port configuration
SERIAL_PORT = '/dev/ttyACM0'  # Replace with your serial port
BAUD_RATE = 9600

# TCP server configuration
TCP_HOST = '0.0.0.0'  # Replace with your server IP if not localhost
TCP_PORT = 6002
datu = "31.2861824,75.5935974"
# Function to read from serial and send over TCP
def serial_to_tcp():
    try:
        # Open serial port
        ser = serial.Serial(SERIAL_PORT, BAUD_RATE)
        print(f"Serial port {SERIAL_PORT} opened successfully.")

        # Create TCP server socket
        tcp_server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        tcp_server_socket.bind((TCP_HOST, TCP_PORT))
        tcp_server_socket.listen(1)
        print(f"TCP server listening on {TCP_HOST}:{TCP_PORT}")

        # Accept incoming connection
        conn, addr = tcp_server_socket.accept()
        print(f"Connection established from {addr}")

        # Function to handle sending data from serial to TCP
        def send_serial_data():
            global datu
            while True:
                data = ser.readline()
                if data.strip() != b'':

                    conn.sendall(data)
                    print(f"Sent data over TCP: {data.strip()}")
                else:
                    conn.sendall(b'31.2861824,75.5935974')


        # Start a separate thread for sending serial data
        serial_thread = threading.Thread(target=send_serial_data, daemon=True)
        serial_thread.start()

        # Main thread continues to handle incoming data or can perform other tasks

        while True:
            pass  # Add your main server logic if needed

    except serial.SerialException as e:
        print(f"Error opening serial port: {e}")
        sys.exit(1)
    except socket.error as e:
        print(f"Socket error: {e}")
        sys.exit(1)
    finally:
        if 'ser' in locals() and ser.is_open:
            ser.close()
        if 'conn' in locals():
            conn.close()
        if 'tcp_server_socket' in locals():
            tcp_server_socket.close()

if __name__ == "__main__":
    serial_to_tcp()
