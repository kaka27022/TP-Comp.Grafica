# TP1 – Visualização 2D de Árvores Arteriais  
**BCC327 – Computação Gráfica**

Este projeto implementa um visualizador interativo de modelos arteriais 2D armazenados em arquivos **VTK (formato Legacy ASCII)**.  
O sistema permite visualizar a árvore, aplicar transformações geométricas, acompanhar o crescimento incremental e ativar recorte (clipping) de segmentos.

---

##  Funcionalidades

###  Leitura de arquivos VTK (Legacy ASCII)
- Suporte a:
  - `POINTS`
  - `LINES`
  - `CELL_DATA` / `POINT_DATA` (raio opcional)
- Reconstrução da geometria a partir de pares de pontos.

###  Transformações geométricas (via matrizes homogêneas 3×3)
- **Translação:** ← → ↑ ↓  
- **Rotação:** `r` (horário), `R` (anti-horário)  
- **Escala:** `+` aumenta, `-` diminui  

###  Crescimento incremental da árvore
- Slider permite exibir os segmentos progressivamente.  
- Também pode ser controlado por:
  - `PageUp` → adiciona 1 segmento  
  - `PageDown` → remove 1 segmento  

###  Recorte (clipping) com algoritmo **Cohen–Sutherland**
- Ativado/desativado pela tecla `c`
- Exibe retângulo de clipping no canvas
- Apenas segmentos dentro da região são renderizados

###  Projeção ortográfica 2D
- Modelos VTK são carregados como coordenadas 2D (descartando Z).

---

##  Interface

Ao executar o programa, você verá:

- A árvore arterial 2D
- Um slider de crescimento na parte inferior
- Um botão “Reset”
- Indicadores de estado:
  - translação atual
  - rotação (em graus)
  - escala
  - número de segmentos exibidos
  - clipping ON/OFF

---

##  Como executar

### 1) Criar ambiente virtual (recomendado)

```bash
python3 -m venv venv
source venv/bin/activate   # Linux/macOS
# venv\Scripts\activate    # Windows
```

### 2) Instalar dependências

```bash
pip install numpy matplotlib
```

### 3) Executar o programa

```bash
python3 tp1_visualizador2d.py
```

---

## Controles

```markdown
| Ação                       | Tecla           |
|----------------------------|-----------------|
| Translação X               | ← →             |
| Translação Y               | ↑ ↓             |
| Rotacionar horário         | r               |
| Rotacionar anti-horário    | R               |
| Aumentar escala            | +               |
| Diminuir escala            | -               |
| Ativar/desativar clipping  | c               |
| Incrementar segmentos      | PageUp          |
| Decrementar segmentos      | PageDown        |
| Reset geral                | Botão “Reset”   |
```

---

## Implementação

### Transformações

As transformações são implementadas usando matrizes homogêneas 3×3:

* Translação

* Rotação

* Escala

Aplicadas na ordem:

```bash
P' = T · R · S · P
```

### Recorte (clipping)

Foi utilizado o algoritmo clássico Cohen–Sutherland, com:

* códigos de região

* testes triviais de aceitação e rejeição

* interseção com bordas