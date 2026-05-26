#!/usr/bin/env python3
import http.server
import socketserver
import webbrowser
import threading
import time
import sys
import os

# Serve from the project root (one level up from scripts/)
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

def find_free_port(start_port=8000):
    port = start_port
    while port < 9000:
        try:
            with socketserver.TCPServer(("", port), None) as s:
                return port
        except OSError:
            port += 1
    return None

def main():
    port = find_free_port()
    if not port:
        print("Error: Could not find an available port to start the web server.")
        sys.exit(1)
        
    class ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
        pass

    Handler = http.server.SimpleHTTPRequestHandler
    # Serve from the directory where the script is located
    server = ThreadingHTTPServer(("", port), Handler)
    
    # Start server in background thread
    server_thread = threading.Thread(target=server.serve_forever, daemon=True)
    server_thread.start()
    
    url = f"http://localhost:{port}/index.html"
    print("=" * 60)
    print("ALPHAGenome UMAP Explorer Server")
    print("-" * 60)
    print(f"Server started successfully!")
    print(f"URL: {url}")
    print("Press Ctrl+C to stop the server.")
    print("=" * 60)
    
    # Open the browser
    try:
        webbrowser.open(url)
    except Exception as e:
        print(f"Note: Could not automatically open the browser ({e}).")
        print(f"Please open your browser manually and go to: {url}")
        
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\nStopping server...")
        server.shutdown()
        server.server_close()
        print("Server stopped.")

if __name__ == "__main__":
    main()
