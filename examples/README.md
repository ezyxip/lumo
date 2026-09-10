# Примеры двух-уровневой задачи

Формализация из обсуждения: два уровня решения.

## Уровень 1 — компилятор: `ColumnConfig -> Model`

Вход компилятора = `compiler_input.json`. Фиксирует **структуру** задачи, меняет размерность/связность MESH, требует перегенерации модели.

Поля:
- `components` — список компонентов с `M` и константами Антуана `A,B,C` (`lg P = A - B/(C+t)`, P в мм рт.ст., t в °C)
- `thermo.model` — `antoine_raoult` (ideal)
- `column.N` — число теоретических ступеней (включая конденсатор и ребойлер)
- `column.condenser/reboiler` — тип
- `column.feed_stages`, `side_draws` — топология
- `column.pressure_MPa` — давление (const по высоте, зашито как compile-time)

## Уровень 2 — модель: `OperatingParams -> Results`

Вход модели = `model_input.json` (дешёвый вызов без перекомпиляции):
- `feed` — `F`, `z` (mass или mole), `e` (доля отгона), `q=1-e`, `tF`
- `specs` — спеки продукта `yD`, `xW`, `R`, `D` и т.п. — замыкают степени свободы
- `initial_guess` — опционально

Выход модели = `model_output.json`:
- `profiles` — `T, x, y, L, V` по ступеням
- `products` — `D, W, Qd, Qw`
- `diagnostics` — `converged, residual, iterations`

## Конверсия mass <-> mole

Методичка (стр. 6, ф. 4-5):

```
x' = (x/M1) / (x/M1 + (1-x)/M2)          # mass -> mole
x  = (M1*x') / (M1*x' + M2*(1-x'))       # mole -> mass
```

Аналогично для y. `M1, M2` — молярные массы.

## Источник данных

- `01` — Губкин, Прокофьева и др. 2016, `doc/theory/Tekhnol_rasch_rektif_kolon_binar_smes.pdf` стр. 4-10, табл. 1-3.
- `02` — бензол/толуол, канонический идеальный пример для кросс-чека с Raoult.
- `03` — синтетический с широкой разницей температур кипения для теста обусловленности.

Все примеры — **бинарные, Antoine+Raoult, P=const**, чтобы reference-солвер считал их правильно.
