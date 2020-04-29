# Лабораторная работа № 6. Безусловный экстремум.

Выполнил студент группы 428  
Мунин Сергей Александрович

## Вариант № 16
Найти точку **максимума**

![FUNC](FUNC.png)

![Pribl](Pribl.png)

методом **наискорейшего спуска**. Для одномерной минимизации использовать метод **квадратичной интерполяции**.
Для поиска интервала унимодальности использовать алгоритм **Свенна**.
В окрестности точки максимума построить линии уровня и траекторию поиска (на одном графике).
## Теоретическая часть

Задаём начальное приближение ![](Pribl.png)
Находим формулы компонент градиента функции в произвольной точке<br>
Далее с помощью алгоритма Свенна, подав на вход приближение ![](Pribl.png) и выбрав шаг, строим широкий интервал. содержащий точку экстремума<br>
### Алгоритм Свенна
В данном методе используется Используется эвристический подход в котором
x
k+1 пробная точка определяется по рекуррентной формуле
![](1.png)
где
![](2.png) - произвольно выбранная начальная точка;
h - шаг поиска, знак которого может меняться на противоположный.
Знак h определяется путем сравнения значений 
![](3.png)
 ![](4.png).<br>
 Если ![](5.png),то согласно предположению об унимодальности, точка минимума должна располагаться правее точки ![](2.png) и величина h выбирается положительной.
Если ![](6.png),то величину h следует выбирать отрицательной. Если ![](7.png),то точка
минимума лежит между ![](8.png) и поиск граничных точек
завершен в противном случае изменить начальную точку. Случай,
когда ![](9.png), противоречит предположению об унимодальности. Выполнение этого условия говорит о том,
что функция в орестности точки ![](2.png) не является унимодальной и
следует изменить начальную точку 

### Метод квадратичной интерполяции

Здесь задаются пробные три пробные точки , одна из которых высчитывается по формуле ![](10.png), как серединная точка интервала из алгоритма Свенна, ![](11.png) и
![](12.png). Для нахождения точки![](11.png)
задается шаг h > 0 в положительном
направлении от точки![](13.png)
, т.е. ![](14.png) и если
![](15.png), то![](16.png), иначе ![](17.png).<br>
Вычисляются значения функции в этих точках ![](18.png),
строится квадратичный интерполяционый многочлен по трем точкам и находится его точка минимума по формуле<br>
![](formula.png)<br>.
Если знаменатель в формуле для нахождения минимума квадратичного интерполяционного многочлена равен нулю, т.е. все три
точки лежат на одной прямой рекомендуется выбрать за ![](20.png)
и повторить нахождение точки ![](21.png)
.
Критерием окончания в этого процесса является выполнение условий для заданного ϵ
![](22.png)
Если условия окончания не выполняются и
![](23.png)
точка ![](13.png)
заменяется на точку ![](24.png), в противном случае точка ![](13.png)
заменяется ![](21.png)

### Метод наискорейшего спуска
 Метод наискорейшего спуска. В этом методе αk выбирается из условия минимума функции ![](25.png) вдоль направления ![](26.png)<br>
, т.е.
![](27.png).<br>
Таким образом, в методе наискорейшего спуска на каждом шаге
необходимо решать задачу минимизации функции одной переменной, что и выполняется в ранее описанных методах.

## Практическая часть
Моя работа состоит из 2-х программ: <br>
1) Программа с реализацией поиска экстремума ф-ции, написанная на языке C++<br>
2) Программа для отрисовки линий уровня ф-ции и траектории поиска экстремума, написанная на Python<br><br>

Программа на C++ состоит из 1-го файла **`Labb6.cpp`**<br><br>

Структура программы:<br>
* В начале программы продключаются библиотеки: <br>
     `iostream` - стандартная библиотека ввода/вывода<br>
     `cmath` - стандартная библиотека для выполнения математических операций <br>
     `fstream` - библиотека для чтения/записи данных из/в файл <br>
* `double F(double x1, double x2)` - функция 2-х переменных, которая дана в билете<br>
* `double grad1(double x1, double x2)` - 1-ая компонента вектора градиента исходной ф-ции<br>
* `double grad2(double x1, double x2)` - 2-ая компоненту вектора градиента исходной ф-ции<br>
* В классе `Sopr` описаны переменные и методы, необходмые для отыскания точки экстремума:
  * *private* - методы и переменные класса:<br>
    * `h` - произвольный шаг, используемый в методе квадратичной интерполяции и методе нахождения интервала унимодальности
    * `k` - счетчик итераций в методе сопряженных градиентов
    * двумерные статические массивы `x0` и  `x1` для записи предыдущего и последующего приближения точки экстремума функции
    * двумерные статические массивы `p0` и  `p1` для записи предыдущего и последующего значений вектора направления                            ![FUNC](Images/vectP.png) в методе сопряж-х град-в
    * `fout` - объект класса ofstream для записи в файл новых значений вектора  ![](Images/xK.png) 
    * ` double f(double alfa)` - та же функция `double F(double x1, double x2)`, только для случая, когда переменные x1 и x2 в свою очередь зависят от  ![](Images/alfa.png). Она необходима для поиска экстремума ф-ции одной переменной в методе квадратичной интерп-ции
    * `double interpMnogochlen(double& x1, double& x2, double& x3)` - квадратичный интерполяционный многочлен, используемый в методе квадр-й интерполяции
    * 'double argmin(double x1, double x2, double x3)' - возвращает один из 3-х аргументов, при котором исходная ф-ция имеет наименьшее значение
    * 'double argmin(double x1, double x2)' - возвращающает один из 2-х аргументов, при котором исходная ф-ция имеет наименьшее значение
    
  * *public* - методы и переменные класса:
    * `int GetIter()` - возвращающает кол-во итераций, за которое удалось найти минимум нашей функции
    * `double* soprGrad()` - реализует *метод сопряженных градиентов*. Возвращает адрес найденной точки минимума.
    * `double KvadrInterp(double a, double b)` - реализует *метод квадратичной интерполяции*. Возвращает минимум ф-ции одной переменной.
    * `void unimodal(double x0, double& a, double& b)` - находит интервал унимодальности нашей ф-ции
   
* В методе `int main()` создается объект типа `Sopr`, у которого вызывается ф-я `soprGrad()`. В консоль выводится найденная точка минимума и кол-во итераций.  
Также ф-я `soprGrad()` создает в папке с программой файл `tr.txt`, в который записываются значения вектора ![](Formuls/xK.png) для построения траектории поиска.<br>
И в конце с помощью команды `system("linur.py")` вызывается питоновский файл, в котором построенна траектория поиска.<br><br>

Программа на Python состоит из 1-го файла **`linur.py`**<br>

В ней строится картина линий уровня нашей ф-ции с помощью метода `contour()` из библиотеки `pylab` и траектория поиска по данным из файла `tr.txt`.<br>

**Порядок компилляции/запуска:**<br>
1. Компиллируем и запускаем файл **`cod.cpp`** из командной строки при помощи команды:<br>
`g++ cod.cpp -o cod.o`<br><br>
2. Запускаем из командной строки `cod.o`<br><br>

### Результаты
В результате работы программы у функции <br><br> ![](FUNC.png) <br><br>был найден экстремум в точке  ![](-6.5,-2) (начальная точка  ![](Images/nachTochk.png)) за ***1*** итерацию методом наискорейшего спуска
 ![](lines.jpg) 
