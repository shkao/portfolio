#!/bin/bash
set -euo pipefail


curl --location --fail https://service.azul.data.humancellatlas.org/manifest/files/ksQwlKVkY3AzOaRjdXJsxBAQ8SWTo8FTpKjSwk0vFvAdxBA_kI4_i05VQ7H5Zurcs6oqxCBLMmxG5y7P6p0HM3uXlJKjETPuFIthphnhWe851fJL3Q \
  | curl --fail-early --continue-at - --retry 15 --retry-delay 10 --config -
