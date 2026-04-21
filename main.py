"""
Точка входа: запускает FastAPI (бэкенд) и Streamlit (фронтенд).

    python main.py

Бэкенд:  http://127.0.0.1:8000
Фронтенд: http://localhost:8501
"""

import subprocess
import sys
import os
import signal
import time

ROOT = os.path.dirname(os.path.abspath(__file__))

# Используем Python из venv если он есть
_venv_python = os.path.join(ROOT, "venv", "Scripts", "python.exe")  # Windows
if not os.path.exists(_venv_python):
    _venv_python = os.path.join(ROOT, "venv", "bin", "python")       # Linux/macOS
PYTHON = _venv_python if os.path.exists(_venv_python) else sys.executable


def main():
    backend_cmd = [
        PYTHON, "-m", "uvicorn", "main:app",
        "--host", "127.0.0.1",
        "--port", "8453",
        "--reload",
    ]

    frontend_cmd = [
        PYTHON, "-m", "streamlit", "run",
        os.path.join(ROOT, "frontend", "Главная.py"),
        "--server.port", "8501",
    ]

    print("Запуск бэкенда  → http://127.0.0.1:8453")
    print("Запуск фронтенда → http://localhost:8501")
    print("Для остановки нажмите Ctrl+C\n")

    backend = subprocess.Popen(
        backend_cmd,
        cwd=os.path.join(ROOT, "backend"),
    )

    time.sleep(1.5)

    frontend = subprocess.Popen(
        frontend_cmd,
        cwd=ROOT,
    )

    def shutdown(sig, frame):
        print("\nОстановка...")
        frontend.terminate()
        backend.terminate()
        frontend.wait()
        backend.wait()
        sys.exit(0)

    signal.signal(signal.SIGINT, shutdown)
    signal.signal(signal.SIGTERM, shutdown)

    # Ждём завершения любого из процессов
    while True:
        if backend.poll() is not None:
            print("Бэкенд завершился неожиданно.")
            frontend.terminate()
            break
        if frontend.poll() is not None:
            print("Фронтенд завершился неожиданно.")
            backend.terminate()
            break
        time.sleep(1)


if __name__ == "__main__":
    main()
