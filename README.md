# 2D модель Изинга: Монте-Карло (Метрополис)

Проект для конференции студенческих научных работ по физике-2026 от ФФ НИУ ВШЭ.  
Реализован на C++ с визуализацией решётки через HTML.

---


## 🚀 Как запустить

### 1. Склонируйте, скомпилируйте и запустите C++ код

```bash
git clone https://github.com/Helen11-png/Ising-model-on-a-square-lattice-Monte-Carlo-method.git
cd Ising-model-on-a-square-lattice-Monte-Carlo-method
g++ -std=c++11 ising_model.cpp -o ising
./ising
```
### 2. После выполнения появятся файлы (это может занять некоторое время):
- results.txt
- lattice.txt
### 3. Откройте view_lattice.html в любом браузере. Выберите файлы lattice.txt и results.txt. Нажмите на кнопку "Построить"
Примечание: визуализация работает локально, без сервера. Файлы lattice.txt и results.txt должен лежать в той же папке, что и view_lattice.html.
