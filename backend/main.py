from fastapi import FastAPI, HTTPException
from schemas.models import FactorizeRequest, FactorizeResponse

from algorithms.pollard import PollardRho
from algorithms.pollard_p1 import PollardP1
from algorithms.quadratic_sieve_basic import QuadraticSieveBasic
from algorithms.quadratic_sieve_optimized import QuadraticSieveOptimized
from algorithms.quadratic_sieve_auto import QuadraticSieveAuto
from algorithms.quadratic_sieve_lpv import QuadraticSieveLPV
from algorithms.quadratic_sieve_mpqs import QuadraticSieveMPQS
from algorithms.quadratic_sieve_mpqs_parallel import QuadraticSieveMPQSParallel

import time
import asyncio

app = FastAPI(title="ВКР: API Факторизации")

ALGO_MAP = {
    "pollard_rho":       PollardRho,
    "pollard_p1":        PollardP1,
    "qs_basic":          QuadraticSieveBasic,
    "qs_optimized":      QuadraticSieveOptimized,
    "qs_auto":           QuadraticSieveAuto,
    "qs_lpv":            QuadraticSieveLPV,
    "qs_mpqs":           QuadraticSieveMPQS,
    "qs_mpqs_parallel":  QuadraticSieveMPQSParallel,
}

TIMEOUT_SECONDS = 30.0

@app.post("/api/factorize", response_model=FactorizeResponse)
async def factorize(request: FactorizeRequest):
    try:
        n = int(request.number)
        if n < 2:
            raise ValueError("Число должно быть больше 1")
    except ValueError:
        raise HTTPException(status_code=400, detail="Некорректный ввод числа")

    algo_cls = ALGO_MAP.get(request.algorithm)
    if algo_cls is None:
        raise HTTPException(status_code=400, detail="Неизвестный алгоритм")

    algo = algo_cls()
    start_time = time.perf_counter()

    try:
        raw_factors = await asyncio.wait_for(
            asyncio.to_thread(algo.factorize, n),
            timeout=TIMEOUT_SECONDS
        )
        factors = [str(f) for f in raw_factors]
        steps = algo.steps_log
        execution_time = (time.perf_counter() - start_time) * 1000
    except asyncio.TimeoutError:
        execution_time = TIMEOUT_SECONDS * 1000
        factors = [str(n)]
        steps = [{"step": "Тайм-аут", "details": {
            "message": f"Превышено время ожидания ({int(TIMEOUT_SECONDS)}с). Алгоритм не справился с числом данной разрядности."
        }}]
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

    return FactorizeResponse(
        factors=factors,
        time_ms=execution_time,
        steps=steps
    )

if __name__ == "__main__":
    import uvicorn
    # Запускаем сервер на порту 8000
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)