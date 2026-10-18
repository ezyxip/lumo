# Общие проверки сгенерированного C99

Пакет используется Go-тестами внутри `codegen`
Он напрямую вызывает системный `cc`, записывает исходники и общий `ctest.h`,
собирает все переданные файлы `.c` и запускает проверки по имени

```go
program := ctest.Build(t, map[string][]byte{
    "model.h": header,
    "model.c": source,
    "check.c": checks,
})
t.Run("balance", func(t *testing.T) {
    program.Run(t, "balance")
})
```

Имена файлов должны быть без вложенных путей
`ctest.h`, `check`, `build.sh` и `build.log` зарезервированы

## Проверки на C

`ctest.h` предоставляет `CTEST_NEAR` для чисел и `CTEST_STATUS` для статусов
Макросы вычисляют аргументы один раз и возвращают 1 из проверки при ошибке
Сообщение содержит имя проверки, файл, строку и сравниваемые значения
Численный допуск равен `1e-10 * max(1, abs(expected))`
Для другого допуска следует добавить явную проверку, а не ослаблять общий

```c
#include "ctest.h"

static int balance(void) {
    CTEST_NEAR(1.0, 1.0);
    CTEST_STATUS(0, 0);
    return 0;
}

int main(int argc, char **argv) {
    const struct ctest_case cases[] = {{"balance", balance}};
    return ctest_run(argc, argv, cases, sizeof cases / sizeof cases[0]);
}
```

`ctest_set_name` позволяет добавить контекст, например температуру
`RunFailure` проверяет ожидаемое падение со статусом 1 и нужными фрагментами
сообщения, он используется в тестах самой инфраструктуры

## Файлы для отладки

Сборка использует `-std=c99 -g -O0 -Wall -Wextra -Werror -pedantic` и `-lm`
При падении сохраняются исходники, бинарник `check`, скрипт `build.sh`,
журнал сборки `build.log` и журналы отдельных проверок
Путь выводится как `C debug artifacts kept in`
После успешного прогона папка удаляется, если не задан `LUMO_CODEGEN_KEEP=1`

```bash
LUMO_CODEGEN_KEEP=1 go test -count=1 -v ./codegen/...
```

Сохранённый `build.sh` повторяет сборку с теми же аргументами
Для GDB можно запустить `gdb --args <путь>/check balance`
После отладки сохранённую папку можно удалить вручную
