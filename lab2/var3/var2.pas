uses crt;

const mf=72;
sigma=5.669e-8;
eps=1e-5;

type
vector1=array[1..mf] of real;
vector2=array[1..mf,1..mf] of real;

var {раздел описания переменных, которые мы будем использовать в
программе}
i, j, Nx, Ny, N1x, N1y: integer;
T, Tn: vector2;
alfa, beta: vector1;
ai, bi, ci, fi, d, d1: real;
a, lamda, ro, c: real;
kapa1, kapa2, Te1, Te2, eps1, q: real;
hx, hy, tau, t_end, time: real;
T0, L, H, alpha: real;
fx, fy, fz: text;

begin
    clrscr;
    
    Nx := 50;
    Ny := 50;
    t_end := 36000;
    L := 0.6;
    H := 0.4;
    lamda := 0.16;
    ro := 1190.0;
    c := 1900.0;
    kapa1 := 50.0;
    kapa2 := 35.0;
    Te1 := 200.0;
    Te2 := 200.0;
    eps1 := 0.8;
    T0 := 300.0;
    q := 10000;
    alpha := 0.1; 
    
    {определяем расчетные шаги сетки по пространственным координатам}
    N1x:=46;
    N1y:=20;
    hx:=L/(Nx-1);
    hy:=H/(Ny-1);
    
    {определяем расчетный шаг сетки по времени}
    tau:=t_end/1000.0;

    {определяем коэффициент температуропроводности}
    a:=lamda/(ro*c);

    for i:= 1 to Nx do
        for j:= 1 to Ny do
            T[i,j]:=T0;
    
    {проводим интегрирование нестационарного уравнения теплопроводности}
    time:=0;
    while time<t_end do
    begin
        time:=time+tau;
        
        {запоминаем поле температуры на n-ом временном слое}
        for i:=1 to Nx do
            for j:=1 to Ny do
                Tn[i,j]:=T[i,j];
        
        {решаем СЛАУ в направлении оси Ох для определения поля температуры на промежуточном (n+1/2) временном слое}
        for j:=1 to Ny do
        begin
            {определяем alfa начальный прогоночный коэффициент на основе левого 
            граничного условия, используя соотношение (48)}
            alfa[1]:=2.0*tau*lamda/(2.0*tau*(lamda+kapa1*hx)+ro*c*sqr(hx));

            {цикл с постусловием, позволяющий итерационно вычислять поле
            температуры, вследствие наличия нелинейности в левом граничном
            условии}
            repeat
                {определяем beta начальный прогоночный коэффициент на основе
                левого граничного условия, используя соотношение (48), при этом
                начинаем итерационный цикл по левому граничному условию}
                d:=T[1,j];
                beta[1]:=(ro*c*sqr(hx)*Tn[1,j]+2.0*tau*kapa1*hx*Te1+2.0*tau
                    *eps1*sigma*hx*(sqr(sqr(Te1))-sqr(sqr(d))))
                    /(2.0*tau*(lamda+kapa1*hx)+ro*c*sqr(hx));
                
                {цикл с параметром для определения прогоночных коэффициентов по
                формуле (8)}
                for i:= 2 to Nx-1 do
                begin
                    {ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с
                    трехдиагональной матрицей}
                    ai:=lamda/sqr(hx);
                    bi:=2.0*lamda/sqr(hx)+ro*c/tau;
                    ci:=lamda/sqr(hx);
                    fi:=-ro*c*Tn[i,j]/tau;

                    if (i=N1x) then fi:= fi - q * Exp(-alpha * (sqr(i - N1x)));
                    //fi:=-ro*c*Tn[i,j]/tau-hx*(i-1)*q;

                    {alfa[i], beta[i] – прогоночные коэффициенты}
                    alfa[i]:=ai/(bi-ci*alfa[i-1]);
                    beta[i]:=(ci*beta[i-1]-fi)/(bi-ci*alfa[i-1]);
                end;
                
                {определяем значение температуры на правой границе, используя
                соотношение (21) при условии, что q 2 = 0}
                T[Nx,j]:=(ro*c*sqr(hx)*Tn[Nx,j]+2.0*tau*lamda*beta[Nx-1])
                        /(ro*c*sqr(hx)+2.0*tau*lamda*(1-alfa[Nx-1]));
                //T[Nx,j]:=(lamda*sqr(hx)*T[Nx,j]+2.0*a*tau*(lamda*beta[Nx-1]+kapa1*hx*Te1)) /(lamda*sqr(hx)+2.0*a*tau*(lamda*(1-alfa[Nx-1])+kapa1*hx));

                {используя соотношение (7) определяем неизвестное поле температуры
                на промежуточном (n+1/2) временном слое}
                for i:= Nx-1 downto 1 do
                    T[i,j]:=alfa[i]*T[i+1,j]+beta[i];
                
            until abs((d-T[1,j]) / T[1, j])<=eps; {значение температуры на левой границе определили}
        
        end; {поле температуры на промежуточном (n+1/2) временном слое определили}
    
        {решаем СЛАУ в направлении оси Оу для определения поля температуры на целом (n+1) временном слое}
        for i:=1 to Nx do
        begin
        {определяем начальные прогоночные коэффициенты на основе нижнего
        граничного условия, используя соотношения (20) при условии, что
        q 1 = 0}
            alfa[1]:=2.0*tau*lamda/(2.0*tau*lamda+ro*c*sqr(hy));
            beta[1]:=ro*c*sqr(hy)*T[i,1]/(2.0*tau*lamda+ro*c*sqr(hy));

            //alfa[1]:=2.0*a*tau*lamda/(lamda*sqr(hy)+2.0*a*tau*(lamda+kapa1*hy));
            //beta[1]:=(lamda*sqr(hy)*T[i,1]+2.0*a*tau*kapa1*hy*Te1)/(lamda*sqr(hy)+2.0*a*tau*(lamda+kapa1*hy));
            
            {цикл с параметром для определения прогоночных коэффициентов по
            формуле (8)}
            for j:= 2 to Ny-1 do
            begin
                {ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с
                трехдиагональной матрицей}
                ai:=lamda/sqr(hy);
                bi:=2.0*lamda/sqr(hy)+ro*c/tau;
                ci:=lamda/sqr(hy);
                fi:=-ro*c*T[i,j]/tau;
                
                //if (j=N1y) then fi:=-ro*c*T[i,j]/tau-hy*(j-1)*q;
                if (j=N1y) then fi:= fi - q * Exp(-alpha * (sqr(j - N1y)));

                {alfa[j], beta[j] – прогоночные коэффициенты}
                alfa[j]:=ai/(bi-ci*alfa[j-1]);
                beta[j]:=(ci*beta[j-1]-fi)/(bi-ci*alfa[j-1]);
            end;
            
            {запоминаем значение температуры на правой границе с
            промежуточного (n+1/2) временного слоя}
            d:=T[i,Ny];
            {цикл, позволяющий итерационно вычислить значение
            температуры на правой границе, вследствие наличия нелинейности в
            этом граничном условии}
            repeat
                d1:=T[i,Ny];
                {определяем значение температуры на правой границе на основе правого граничного условия}
                T[i,Ny]:=(ro*c*sqr(hy)*d+2.0*tau*(lamda*beta[Ny-1]+kapa2*hy*Te2
                +eps1*sigma*hy*(sqr(sqr(Te2))-sqr(sqr(d1)))))/(ro*c*sqr(hy)
                +2.0*tau*(lamda*(1-alfa[Ny-1])+kapa2*hy));
            until abs((d1-T[i,Ny]) / T[i, Ny])<=eps; 
            
            for j:= Ny-1 downto 1 do
                T[i,j]:=alfa[j]*T[i,j+1]+beta[j];
        end; {поле температуры на целом (n+1) временном слое определили}
    end;

    Writeln('Длина пластины L = ', L:0:2);
    Writeln('Ширина пластины H = ', H:0:2);
    Writeln('');
    Writeln('Число узлов по x = ', Nx);
    Writeln('Число узлов по y = ', Ny);
    Writeln('');
    Writeln('Начальная температура T0 = ', T0:0:3);
    Writeln('Температура окружающей среды T1 = ', Te1:0:4);
    Writeln('');
    Writeln('Шаг по x = ', hx:0:3);
    Writeln('Шаг по y = ', hy:0:3);
    Writeln('Время t = ', t_end:0:3);

    Assign(fx, 'X.txt');
    Rewrite(fx);
    for i:=1 to Nx do
        Writeln(fx, hx * (i-1):0:8);
    close(fx);

    Assign(fy, 'Y.txt');
    Rewrite(fy);
    for j:=1 to Ny do
        Writeln(fy, hy * (j-1):0:8);
    close(fy);

    Assign(fz, 'Z.txt');
    Rewrite(fz);
    for i:=1 to Nx do
    begin
        for j:=1 to Ny do
            Write(fz, T[i,j]:0:8, ' ');
        Writeln(fz, '')
    end;
    close(fz);

end.
