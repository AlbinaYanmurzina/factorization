#!/usr/bin/env bash
# Запуск бэкенда и фронтенда в фоне.
# Использование: ./start.sh
# Остановка:     ./start.sh stop

set -e

ROOT="$(cd "$(dirname "$0")" && pwd)"

# Активируем venv если есть
if [ -f "$ROOT/venv/bin/activate" ]; then
    source "$ROOT/venv/bin/activate"
fi

PIDFILE="$ROOT/.pids"

stop() {
    if [ -f "$PIDFILE" ]; then
        echo "Остановка процессов..."
        while IFS= read -r pid; do
            kill "$pid" 2>/dev/null && echo "  убит PID $pid" || true
        done < "$PIDFILE"
        rm -f "$PIDFILE"
        echo "Готово."
    else
        echo "Нет запущенных процессов (файл .pids не найден)."
    fi
    exit 0
}

if [ "${1:-}" = "stop" ]; then
    stop
fi

# Проверяем, не запущено ли уже
if [ -f "$PIDFILE" ]; then
    echo "Уже запущено (найден .pids). Сначала выполните: ./start.sh stop"
    exit 1
fi

echo "Запуск бэкенда  → http://127.0.0.1:8000"
cd "$ROOT/backend"
uvicorn main:app --host 127.0.0.1 --port 8000 --reload &
BACKEND_PID=$!

sleep 1.5

echo "Запуск фронтенда → http://localhost:8501"
cd "$ROOT"
python -m streamlit run frontend/Главная.py --server.port 8501 &
FRONTEND_PID=$!

# Сохраняем PID-ы для команды stop
printf "%s\n%s\n" "$BACKEND_PID" "$FRONTEND_PID" > "$PIDFILE"

echo ""
echo "Оба процесса запущены."
echo "Для остановки: ./start.sh stop"
echo "  Бэкенд  PID: $BACKEND_PID"
echo "  Фронтенд PID: $FRONTEND_PID"

# Ждём завершения обоих (Ctrl+C остановит оба)
trap 'stop' INT TERM
wait $BACKEND_PID $FRONTEND_PID
