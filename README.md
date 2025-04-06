# ��������� ������������
## �������� �������
� ������ ������� ���� ������������ ��������� ��������� ���������� ���������� ��������� ������������ (��������� � ��� ����).

### ���� ������
�������� �������� ���������� ��������� ������������, ��������� 3 ���� ����������:
* "�������" - ����������, ��������� ���� ������ ������� �� ���� �������� �����;
* "�����������������" - ����������� ���������� �� ���� ��������� ���������� ����� ������������
* SIMD-����������� (AVX2) � ��������� �������� � �������������� ��������� ���������� ����������.

## ��������� ������������
���������� ����� ��������� ������������ ������������� ��������� ����������:
1) $x_{n+1}$ = $x_n^2 - y_n^2 + x_0$
2) $y_{n+1}$ = $2 * x_n * y_n + y_0$,

��� (x0, y0) - ���������� ��������������� �������.  

���� ������� ����������� ������� �� ���������� �������� (`MAX_ITERATIONS = 256`), ��� ������� ����� ��������� ������ ���������� ������� ������� (`MAX_RADIUS = 10`). ����� ��������, ���, ���� ����� ���������� �� ������ ���������� �� ���������� �������, ��� �� ������, �� ������ ��� ����� ��������� � "�������������".

![����](pic/mandlebrot.png)  

����� �������, ����� ����� ��������� �������, � ������� ����� �� ����� �� ������� ���������� (������ ����) � �������� (�� ������ ����).

## ����������� � ������ ���������
### �����������
� ������ ������� �������� ��� ������ ������:
1) **����������� �����**: ����������� ����������� ���� �������� `800x600`, � ������� ��������������� ��������� ������������. � ����� ������� ���� ������������ FPS. ����� �������������� ����������� � ��������� �������� ������.
 
������� �������:
  
|   �������   |   ��������    |
|------------:|:--------------|
| ?           | ����� �����   |
| ?           | ����� ����    |
| ?           | ����� �����   |
| ?           | ����� ������  |
| =           | �����������   |
| -           |  ��������     |  

2) **�������� �����**: ����������� ���� ������� 80 ������. ���������� �������� �������� �� ���������� �������� (80). � ����������� ����� ���������� ����������� ����� ������� ������ �����, ���� � ��������, ���� � ����� ����������
� ����������� �� ������ ���������, ���������� �������������.

� ��������� ����������� ��� ��������� ������ ���������� ������� ��������:
  - **Naive** - ��������� �������� ������� �������.
  - **Array** - ��������� �� 1 �������� ����� ������� �� 8 ��������.
  - **Simd**  - ��������� �� 1 �������� ����� 8 �������� c �������������� AVX ����������, ��� ��������� �������� ������� �����������.

### ������
��� ����, ����� ������� ������, ����� ������ �������:
```bash
make
```

����� ����� ������� �������� `*.exe` ����� � ����� ���� ������� ��������� ��������� ������.  

��� ����, ����� ��������� ����������� ���������� ����� �������, ����� ����� ���������� �� ������ ������������ (1 - `Naive`, 2 - `Array`, 3 - `Simd`).   

������ ������� ����������� "�������" ����������:
```bash
./mandle.exe 1  
```
  
��� ����, ����� ��������� �������� ���������� ����� ������� �������������� �������� ��������� ������, ������� �������� �� ��, � ��� �� ������ �������� ��������� (0 - `seconds`, 1 - `ticks`).  

������ ������� �������� "�������" ���������� � ����������� � ����� ����������:
```bash
./mandle.exe 1 1
```
### ��������� ������������

���������: 12th Gen Intel(R) Core(TM) i5-1200H, 2500 ���, 12 ����, 16 ���������� �����������.

����������: gcc 13.3.0 & clang 18

�C: Ubuntu WSL2

# ������ ������ 
### ������� ���������� ��������� ������������
��������� ������� 1 ����� �� ����� �������� � ������ ����������, ��������� ������� `__rdtsc()`. ������ ����� �������� ���������� ������, ��������� `__rdtsc()` ������������ � ��������� ������������ ���������� (�� ����������� �� ��������� ���������� �������� ������). �������� ����������� � ��������� ����� ������� ����������, ��������� ������� ������ ���������.   

������ ���� �������� ��� ������ ������ ��������� � ������ �� ������� �����������:
 - **-O0**: ����������� ���� �� ������������. �������� ����� ��������� ������������ ���������������� ���������������� � �������� ���. ��������� ����� ����������� ����� ���������� ���������.
 - **-O3**: ������ ����������� ��� ��������� ����������� ������� �����������.
 - **No flag**: ����� ���� �������� ������ ��� ������� ��������� ��� ���������� ���� ������.
  
### ��������� �����������
��� ������ ����� ��������� ������ ��������� � [�������](data.md). ����� ����� ����������� ��������, ��� ��� ��� �������� ����� �������.  

| Type\Flag       | No flags     |    -O0     |    -O3     |
|----------------:|:------------:|:----------:|:-----------|
|Naive on GCC           | 448836152    | 444215145  | 215564328  |
|Array on GCC          | 596859921   | 593526982 | 107144019  |
|Simd on GCC            | 186629892    | 188514684  | 32862219   |  
|                       |              |             |           |
|Naive on Clang         |  447153636   | 446413405   | 210385565 |
|Array on Clang         |  575038952   | 563470040   | 111039172 |
|Simd on Clang          |  182377834   |  183862717  | 31841715  |


��� ������� ��������� ������������� ������ ���������� � ������� �����������.
![����](pic/gis_base.png) 

����� ��� Clang 
![����](pic/gis_base2.png) 

## ���������� ������������ ������������������ ������ ������
### �������� ����������

1. **������������� ����������� -O3**  
   ��� ������ ���������� ������������� ������������ ��������� ��� ���������� � ������ `-O3`:  
   - ��������� � **2.0-5.8?** �� ��������� � `-O0`  
   - ���������� ������� ������������������ ����������� � SIMD on Clang (������� ������������������ ������������ No flags � 5.8 ���)
2. **������������ ������������������ (-O3)**  
   SIMD on GCC:         32.8 ��� ������ (6.5? ������� Naive)
   Vectorized on GCC:   107.1 ��� ������ (2.0? ������� Naive)
   Naive on GCC:        215.5 ��� ������ (������� ��������)

   SIMD on clag:        31.8 ��� ������ (6.6x ������� Naive)
   Vectorized on Clang: 111.0 ��� ������ (1.9x ������� Naive)
   Naive on Clang:      210.3 ��� ������ (������� �������)
3. **��������� ������������**
   
   | �����       | GCC (-O3) | Clang (-O3) | ������� |
   |-------------|----------|------------|---------|
   | Naive       | 215.5M   | 210.3M     | +2.5%   |
   | Vectorized  | 107.1M   | 111.0M     | -3.5%   |
   | SIMD        | 32.8M    | 31.8M      | +3.1%   |

   Clang ���������� marginally ������ ���������� ��� SIMD (+3.1%), GCC ����� ������������ ��������� ������.




### ������ ��� -O3 �������� ������?
������������� **Godbolt** � ���������, �� ��� ������������ ������� c ������������ � ������ -O3 � ��� ����. �� �����, ��� ������ (������ ��� -O3) ������������ � ������� ���������� ����������, � ������ ������������� `vmovaps`, ��� ��� ������������� ������� ���������� ���������. � ������ ����� (������ � -O3) ��������. ������� �� � �������� ������� ���������� ������ ���������� � ������ -O3.
![����](pic/godbolt.png)

# ������ ������
## ����������� ������������������ ����� ���������� ������� �������

���������� � ������� ������������� ������ ������� **8 ���������**, ��� �������������:
- 256-������ ��������� AVX2 (����� 8 float-��������)
- ������������� ������������ ��� ����������� SIMD-����������

**���� ������ � ������������� ��� ���, ����� ���������� ������ ������������**

**��������**: ���������� (GCC/Clang) �� �������� ����������� �����������, ������ ���:
```
float X0[8]; // ���������� ������������ ��� ��������� ������
float X1[8]; // ��� ������� ������������ ������������
```

����������� � ```VECTOR_SIZE = 32```

������ ���������� �������� � ���������, ������� ����������� �������� �������� ������������ �������

### ```On GCC```


  | Type\Flag       |    -O0       |    -O3     |
|----------------:|:------------:|:----------:|
|Array on GCC         | 780228580 | 39062525  |
|Simd on GCC          | 169272774  | 33721441 | 

> [!NOTE]
> ����� ����� ���� ����������� ��� ���������� ������������� ���, �� ����� ��������������� ������ ```-fopt-info-vec```

### ```On Clang```

| Type\Flag       |    -O0       |    -O3     |
|----------------:|:------------:|:----------:|
|Array on Clang        | 851709858 | 91598302  |
|Simd on Clang         | 161623396  | 34310658 | 

## ������������� ��������� 
**GCC** 
  - Array ��� ����� O0 ����������������� ��������� � 0.7 ���, �� ��� ����� O3 �������� ������� � 2.7 ���
  - Simd ��� ����� O0 ������������������ ���������� � 3.5 ���, �� ��� ����� O3 �������� ��������� � 0.97 ���

**Clang**
  - Array ��� ����� O0 ������������������ ���������� � 0.7 ���, �� ��� ����� O3 �������� ������� � 1.2 ����
  - Simd ��� ����� O0 ������������������ ���������� � 1.2 ����, �� ��� ����� O3 �������� ��������� � 0.95 ���

** ���������� **

��� ���������� ������������ � ������ �������

## ������ ��� ���������� 

**��� Array (mandelbrot_vectorized):**

```O0```

- ���������:
    - ���������� �� ������������ ���
    - ������������� ����� ������ � ������ ����� (� 4 ���� ������ ���������)
    - �������� �� ���-������ ����������, �� ��� ����������� ��� �� ��������������

 ```O3```

- ��������� ���������
    - ������� ������������� ���� (������� ��������� ����������)
    - ���������� ��������� �������� �� ����� (������ ������ �� ��������)
 
**��� SIMD(mandelbrot_simd):**

```O0```

- ��������� ������ ���:
  - ���� ��� �����������, ��������� 32 ��������� ������ 8 ����� ���������� ���������� �����������
  - ����������� ���������� �������� �������� �����
 
```O3```

- �������������� ��������� �����������
    - ����������� �������������� ������ CPU (���� saturation)
    - ���������� � ������������� ������ ��� ������� VECTOR_SIZE
    - ����������� �������� �� ���-������

### ����� �� ���� ������

1. ��� ������� ����������:
   - VECTOR_SIZE=32 ������� ������ � -O3 (?2.7x � GCC)
   
2.��� ������ SIMD-�����������:
   - ���������� VECTOR_SIZE=8 (������ IPC)
   - ���������� ������� ���� marginal ������� ������ ��� -O0
 
# ������ ������ 
###  �������� ��������� � ����������

```
 // 1. ����������� ������� ��� ����������
struct MandelBrot {
    alignas(32) float X0[32];
    alignas(32) float X_N[32];
    alignas(32) float Y_N[32];
    // ... ��������� �������
};

// 2. ���������������� SIMD-���������
void mandelbrot_simd(MandelBrot_t* set) {
    __m256i* pixels = (__m256i*)set->pixels_array;
    for (...) {
        __m256 X0 = _mm256_load_ps(&set->X0[0]); // Aligned load
        // ... ���������� ...
        _mm256_store_ps(&output[0], result); // Aligned store
    }
}
```

��� ��������� � ������� ����� ��� � ����������� [```mandelbrot_alignment.cpp```](https://github.com/MaxGud10/Mandelbrot/blob/main/src/mandel_alig.cpp)

** ��������� ������������ **

| ��� ���������� | ���������� | -O0 (�����) | -O3 (�����) | ��������� -O3/-O0 |
|----------------|------------|-------------|-------------|-------------------|
| Array          | GCC        | 541945398 | 105964364 | 5.11x             |
| SIMD           | GCC        | 166633068 | 32932176  | 5.06x             |
| Array          | Clang      | 559299600 | 108979561 | 5.13x             |
| SIMD           | Clang      | 158135311 | 30832159  | 5.13x             |

## ������������� ������

### 1. ������������� ����������� (-O3)

**�������� ����������:**
- ��� ���������� ���������� **��������� ~5.1x** ��� ��������� ����������� (-O3)
- ����������� ������������� ������ �������������:
  - GCC: 5.06-5.11x
  - Clang: 5.13x

**�������������:**

���������������� �������� ��������� ������ � ����� ������������

### 2. ��������� ����������

**Array vs SIMD:**
| �������          | Array (-O3) | SIMD (-O3) | ������������ SIMD |
|------------------|-------------|------------|-------------------|
| GCC              | 105.9M      | 32.9M      | 3.22x             |
| Clang            | 109.0M      | 30.8M      | 3.54x             |

**������:**
- ������ SIMD-����������� ���� **3.2-3.5x ������������** ��� ��������� �������
- ������� ����� ������������� ������������� ( < 10%)

### 3. �������� ��������� ������������ ������� ������

**GCC:**
- **Array**:
  - -O0: +10% (1.1x)
  - -O3: ��� ��������� (1.0x)
  
- **SIMD**:
  - -O0: +10% (1.1x)
  - -O3: -1% (0.99x)

**Clang:**
- **Array**:
  - -O0: ��� ��������� (1.0x)
  - -O3: ��� ��������� (1.0x)
  
- **SIMD**:
  - -O0: +20% (1.2x)
  - -O3: ��� ��������� (1.0x)

# ��������� ������
## ������������� ����������

### ��������� �������

**��������������** ����������
- **����**: 12 ����������  
- **������**: 16 

**���������** ������� ������������������
  - **���**: 12-������� ��������� �� ���� 12 ����  
  - **����������**: ������� ������� �� ���� ����� � ���������� ���������������  

����������� ���������� ������������ ���������� ��������������� (Hyper-Threading � Intel ��� SMT � AMD), ������� ��������� ������ ����������� ���� ��������� ����� ��� ����������� ������.

�������� ���� ��������������� ����������� � ����� ����������� ������������� �������������� ��������: ���� ���� ����� ��� ������ �� ������ ��� �� �����-�� ���� ������� �����������, ������ ����� ����� ������������ ��������� �������������� ����� ����

��� ������ �������� ���������� ��������������� ����� ���������� �������������� ��������� ������� 20�30% �� ��������� � �������������� ������ ���������� ���� (x7-8 � ����� ������). �� ����� ���������� � ������� � � ����� ������, ����� ����� ������ �������� �������������� ��������, ������� ����� ��������� ����, ��� ��� �����, ��� ����� ��������� ��������.

## ������������� OpenMD

> [!NOTE]
> OpenMP (Open Multi-Processing) � �������� �������� ��� ����������������� ��������. ��� �������� ������������ �������� ����������� (#pragma ...), ������������ �������� � ���������� ���������, ������� ������������� ��� ���������������� ������������� ���������� �� ����������������� �������� � ����� �������.

## ���������� ������ �����������������
**OpenMP** ���������� ��������������� ������ - ���������� ����� ���������:

```
    #pragma omp parallel for
    for (size_t y = 0; y < HEIGHT; y++) 
    {
        const float base_x0 = (-(HALF_WIDTH) * D_X + set->x_offset) * set->scale;
    . . . 
```

## ���������� ������ �������

**������������� ��������** - �������� ����������:
  - ����� c_y ����������� ����� ��������� (����������� ����� ����������)
  - ����� ������ ��������� ���������� ��������� ��� ������� Y

## �������� ���� ������������� ��������

### 1. `schedule(dynamic)`
- **������� ������**:  
  - ������������ ������������� �������� �� ���� ���������� �������
  - ������ ����� (����� ��������) �� ��������� = 1
- **������������**:
   - ������� ��� ������������� ��������
   - �������������� ������������  
- **����������**:
   - �������������� ��������� �������  

### 2. `schedule(dynamic, N)`
- **�����������**:  
  - ������������� ������ ����� = N ��������
  - ������ �������� ����� �� 8 ��������
- **����������� ����������**:
    - ����� ����� ���������� �������� �������� ���������  
    - ��������� ��������� ������� vs `dynamic`  

### 3. `schedule(guided, 8)`
- **��������� ��������**:  
  - �������� � ������� ������ ? ��������� ������ �� ���� ����������
  - ����������� ������ ����� = 8 (��������� �����)
- **������������**:  
  - ������� ����� (������� �����)  
  - ������� ������������ (������ ����� � �����)  

## ���������� 

| Type\Flag       | Array        |   SIMD     | 
|----------------:|:------------:|:----------:|
|Dynamic          | 35954467    | 9213873  |
|Dynamic, N        | 38602703   | 9269230 | 
|Guided           | 35969664    | 8758023  |  

## �������� ����������

### 1. ��������� ������� ������������
- **������ ���������**:  
  - ��� Array: **Static** ������� � 3x ������������ ������ ������  
  - ��� SIMD: **Guided** ������� � 3.7� ������������ ������ ������  



# ���������� 

## ������� ������ ��� ���������� � ��� ����������������� ��������

|����������|���������, Ticks| ������� ��� ����� ������ ��� �����������|
|----------------:|:---------------:|:---------------:|
|Naive                             | 210385565    | � 2.0x|
|Array - ������                   |   107144019   | � 5.7x|
|Simd - ������                    | 31841715   | � 5.8x |
|                                 |             | 
|Array - ����������� �� ��������  | 39062525    | � 15.0x |
|Simd - ����������� �� ��������   | 33721441    | � 5.5x |
|                                 |             | 
|Array � �������������            |  105964364  | � 5.6x |
|Simd � �������������             |   30832159  | � 6.0x   |
|                                 |             |
|������������� - dynamic          | 35954467  | �  16.6x|
|������������� - guided           |  8758023  | �  21.3 |