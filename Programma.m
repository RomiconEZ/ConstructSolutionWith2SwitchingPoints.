
clc;
clear all;
close all;
%format long
accuracy = input('Enter accuracy ...');
options = odeset('RelTol',accuracy);
% str ='Y= AY + BF(\sigma) + Kf(t)';
% disp(str)
str ='1) f(t)=f_{0} + f_{1} * sin(wt + \phi_{1}) + f_{2} * sin(2wt + \phi_{2})';
disp(str)
str ='2) ...another type...';
disp(str)
type_of_function = input('Enter the number corresponding to the type of function...');
if type_of_function==1
    f0 = input('Enter f_0 ...');
    f1 = input('Enter f_1 ...');
    f2 = input('Enter f_2 ...');
    w = input('Enter w ...');
    phi1 = input('Enter \phi_1 ...');
    phi2 = input('Enter \phi_2 ...');
    
%     syms t;
%     fs = input('Введите текст функции: ');
%     f = symfun(fs,t);
    
     
    syms t;
    fs = f0 + f1 * sin(w*t + phi1) + f2 * sin(2*w*t + phi2);
    fs = subs(fs,f0,f0);
    fs = subs(fs,f1,f1);
    fs = subs(fs,f2,f2);
    fs = subs(fs,w,w);
    fs = subs(fs,phi1,phi1);
    fs = subs(fs,phi2,phi2);
    f = symfun(fs,t);
    
    fprintf('f(t)=');
    disp(f);
    %%
    raz = input('Enter dimension of the space ...');% размерность пространтсва

    A = zeros(raz);
    for i=1:1:raz
        for j=1:1:raz
           text =['Enter A[',num2str(i),',',num2str(j),'] ...'];
           A(i,j)=input(text);
        end     
    end
    %disp(A);
    B = zeros(raz,1);
    for i=1:1:raz
        text =['Enter B[',num2str(i),',1] ...'];
        B(i,1)=input(text);    
    end
   proverka_uslov_poln_uprav=1;
   for i=0:1:(raz-1)
       if A * B.^i == 0
           proverka_uslov_poln_uprav=0;
       end
   end
   if proverka_uslov_poln_uprav==1
        proverka_obr_preob=1;
        for i=0:1:(raz-1)
            if A^i * B == 0
                proverka_obr_preob=0;
            end
        end
        if proverka_obr_preob==1
            K = zeros(raz,1);
            for i=1:1:raz
                text =['Enter K[',num2str(i),',1] ...'];
                K(i,1)=input(text);    
            end
            C = zeros(raz,1);
            for i=1:1:raz
                text =['Enter C[',num2str(i),',1] ...'];
                C(i,1)=input(text);    
            end
            %T = input('Enter T(period) ...');

            E = eye(raz);
            syms D(p);
            D(p) = A - p*E;
            %disp(D(p));

            syms p;
            AD=cell(raz,raz);% create a cell array
            for i=1:1:raz
                for j=1:1:raz
                   AD{i,j}=ad(i,j,D(p));% ad - алгебраическое дополнение
                end     
            end

            N=cell(raz,1);
            for i=1:1:raz
                N{i,1}(p)=B(1,1)*AD{1,i};
                 for k=2:1:raz
                    N{i,1}(p)=N{i,1}+B(k,1)*AD{k,i};
                 end
            end

            S = solve(det(D(p)),p);%вектор корней характеристического уравнения
            S=sort(S,'descend');
            disp('Корни характеристического уравнения: ');
            disp(S);
            proverka_korney_harak_urav=1;
            for i=1:1:raz
                if S(i,1)==0
                    proverka_korney_harak_urav=0;
                end
                if isreal(S(i,1))==false
                    proverka_korney_harak_urav=0;
                end
            end
            if numel(unique(S))~=raz
                     proverka_korney_harak_urav=0;
            end
            if proverka_korney_harak_urav==1
                syms difD(p);% производная характеристического уравнения
                difD(p) = diff(det(D(p)));
                %disp(difD(p));
                GAMMA=zeros(raz,1);%вектор обратной связи
                for i=1:1:raz
                    sumSN=0;
                    for k=1:1:raz
                        sumSN = sumSN + C(k,1)*double(N{k,1}(S(i,1)));
                    end
                    GAMMA(i,1)=double((-1/difD(S(i,1)))*sumSN);%компоненты вектора обратной связи
                end
                %%
                disp('Г: ');
                disp(GAMMA);
                Smatrix = zeros(raz);%Smatrix -матрица S в статье
                for i=1:1:raz
                    for j=1:1:raz
                       Smatrix(i,j)=double(N{i,1}(S(j,1))/difD(S(j,1)));
                    end     
                end
                Smatrix=Smatrix*(-1);
                disp('Матрица преобразований S: ')
                disp(Smatrix);
                K_0 = Smatrix^(-1) * K;
                disp('K_0: ')
                disp(K_0);
                B_0 = Smatrix^(-1) *B;
                disp('B_0: ');
                disp(B_0);
                A_0 = Smatrix^(-1) *A * Smatrix;
                disp('A_0: ');
                disp(round(A_0,13));%округление до 13 знаков после запятой в силу очень маленьких значений, которые не приводят матрицу к диагональному виду
                gamma_s = GAMMA(1,1);
                i=1;
                while gamma_s==0% находим не нулевую гамму
                    i=i+1;
                    gamma_s = GAMMA(i,1);
                end
                lambda_s=double(S(i,1));
                s_position = i;
                sigma_1 = double(atan(w/lambda_s));
                sigma_2 = double(atan(2*w/lambda_s));

                %H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));
                %syms H(t)
                %H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                %syms L(l1in,k);
                %L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));

                %%
                % myfun = @(x,c) cos(c*x);  % parameterized function
                % c = 2;                    % parameter
                % fun = @(x) myfun(x,c);    % function of x alone
                % x = fzero(fun,0.1)
                %m2=(B_0(s_position,1)*m1*exp(lambda_s*T/k))/(exp(lambda_s*T/k)-exp(lambda_s*t1))+((gamma_s*K_0(s_position,1)*f1*sin(2*pi/k+phi1+sigma_1))/sqrt(lambda_s^2+w^2) + (gamma_s*K_0(s_position,1)*f2*sin(4*pi/k+phi2+sigma_2))/sqrt(lambda_s^2+4*w^2) + l1 + gamma_s*K_0(s_position,1)*f0/lambda_s - exp(lambda_s*T/k)*(l1+gamma_s*(B_0(s_position,1)*m1+K_0(s_position,1)*f0)/lambda_s + gamma_s*K_0(s_position,1)*M))*(lambda_s*exp(lambda_s*t1))/(gamma_s*(exp(lambda_s*T/k)-exp(lambda_s*t1))); 
                fun_po_m2 = @(T,m2,m1,lambda_s,k,gamma_s,f1,f2,f0,phi1,phi2,sigma_1,sigma_2,w,l1,B_0,K_0,M,Ll1k) m2-(B_0(s_position,1)*m1*exp(lambda_s*T/k))/(exp(lambda_s*T/k)-exp(lambda_s*double(T/k + 1/lambda_s*log((m2-m1)/((exp(lambda_s*T/k)-1)*Ll1k+m2-m1*exp(lambda_s*T/k))))))+((gamma_s*K_0(s_position,1)*f1*sin(2*pi/k+phi1+sigma_1))/sqrt(lambda_s^2+w^2) + (gamma_s*K_0(s_position,1)*f2*sin(4*pi/k+phi2+sigma_2))/sqrt(lambda_s^2+4*w^2) + l1 + gamma_s*K_0(s_position,1)*f0/lambda_s - exp(lambda_s*T/k)*(l1+gamma_s*(B_0(s_position,1)*m1+K_0(s_position,1)*f0)/lambda_s + gamma_s*K_0(s_position,1)*M))*(lambda_s*exp(lambda_s*double(T/k + 1/lambda_s*log((m2-m1)/((exp(lambda_s*T/k)-1)*Ll1k+m2-m1*exp(lambda_s*T/k))))))/(gamma_s*(exp(lambda_s*T/k)-exp(lambda_s*double(T/k + 1/lambda_s*log((m2-m1)/((exp(lambda_s*T/k)-1)*Ll1k+m2-m1*exp(lambda_s*T/k))))))); 

                k = 1;
                postr_gr = input('Enter 1 if you want to enter "l1" or 2 if you want to enter "l2" ...');
                    if postr_gr == 1
                        l1 = input('Enter l1 ...');
                        m1 = input('Enter m1 ...');
                        vvod_m2_T = input('Enter 1 if you want to enter "m2" or 2 if you want to enter "T" ...');
                        if vvod_m2_T==1
                            m2 = input('Enter m2 ...');
                            T = input('Enter T(period) ...');
                            H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));       
                            syms H(t)
                            syms L(l1in,k);
                            H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                            L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));
                        else
                            T = input('Enter T(period) ...');
                            m2 = input('Enter m2 ...');
                            %M = (f1*sin(phi1+sigma_1))/sqrt(lambda_s^2+w^2) + (f2*sin(phi2+sigma_2))/sqrt(lambda_s^2+4*w^2);
                            H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));       
                            syms H(t)
                            syms L(l1in,k);
                            H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                            L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));
                            %k=1;
                            %Ll1k = double(L(l1,k)); 
                            %myfun_po_m2 = @(m2) fun_po_m2(T,m2,m1,lambda_s,k,gamma_s,f1,f2,f0,phi1,phi2,sigma_1,sigma_2,w,l1,B_0,K_0,M,Ll1k);
                            %m2 = fsolve(myfun_po_m2,0);
                        end
                        k=1;
                        t1 = double(T/k + 1/lambda_s*log((m2-m1)/((exp(lambda_s*T/k)-1)*L(l1,k)+m2-m1*exp(lambda_s*T/k))));
                        fun = symfun(exp(lambda_s*(t1-t))*f(t),t);
                        fun = subs(fun,lambda_s,lambda_s);
                        fun = subs(fun,t1,t1);
                        intfun = double(int(fun, 0, t1));
                        l2 = double((l1+gamma_s*m1/lambda_s)*exp(lambda_s*t1)-gamma_s*m1/lambda_s+gamma_s*K_0(s_position,1)*intfun);
                    else
                        l2 = input('Enter l2 ...');
                        m1 = input('Enter m1 ...');
                        vvod_m2_T = input('Enter 1 if you want to enter "m2" or 2 if you want to enter "T" ...');
                        if vvod_m2_T==1
                            m2 = input('Enter m2 ...');
                            T = input('Enter T(period) ...');
                            H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));       
                            syms H(t)
                            syms L(l1in,k);
                            H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                            L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));
                        else
                            T = input('Enter T(period) ...');
                            m2 = input('Enter m2 ...');
                            %M = (f1*sin(phi1+sigma_1))/sqrt(lambda_s^2+w^2) + (f2*sin(phi2+sigma_2))/sqrt(lambda_s^2+4*w^2);
                            H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));       
                            syms H(t)
                            syms L(l1in,k);
                            H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                            L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));
                            %k=1;
                            %Ll1k = double(L(l1,k)); 
                            %myfun_po_m2 = @(m2) fun_po_m2(T,m2,m1,lambda_s,k,gamma_s,f1,f2,f0,phi1,phi2,sigma_1,sigma_2,w,l1,B_0,K_0,M,Ll1k);
                            %m2 = fsolve(myfun_po_m2,0);
                        end
                        k=1;
                        t1 = double(T/k + 1/lambda_s*log((m2-m1)/((exp(lambda_s*T/k)-1)*L(l1,k)+m2-m1*exp(lambda_s*T/k))));
                        fun = symfun(exp(lambda_s*(T/k-t))*f(t),t);
                        fun = subs(fun,lambda_s,lambda_s);
                        fun = subs(fun,T,T);
                        fun = subs(fun,k,k);
                        intfun = double(int(fun, t1, T/k));
                        l1 = double((l2+gamma_s*m2/lambda_s)*exp(lambda_s*(T/k-t1))-gamma_s*m2/lambda_s+gamma_s*K_0(s_position,1)*intfun);
                    end
                while (double(m2-m1*exp(lambda_s*T/k)+(exp(lambda_s*T/k)-1)*double(L(l1,k)))<=0) || (m1>=double(L(l1,k))) || (double(L(l1,k))>=m2) || (l1>l2) || (m1>m2)   
                    postr_gr = input('The entered data is incorrect. Enter 1 if you want to enter "l1" or 2 if you want to enter "l2" ...');
                    if postr_gr == 1
                        l1 = input('Enter l1 ...');
                        m1 = input('Enter m1 ...');
                        vvod_m2_T = input('Enter 1 if you want to enter "m2" or 2 if you want to enter "T" ...');
                        if vvod_m2_T==1
                            m2 = input('Enter m2 ...');
                            T = input('Enter T(period) ...');
                            H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));       
                            syms H(t)
                            syms L(l1in,k);
                            H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                            L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));
                        else
                            T = input('Enter T(period) ...');
                            m2 = input('Enter m2 ...');
                            %M = (f1*sin(phi1+sigma_1))/sqrt(lambda_s^2+w^2) + (f2*sin(phi2+sigma_2))/sqrt(lambda_s^2+4*w^2);
                            H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));       
                            syms H(t)
                            syms L(l1in,k);
                            H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                            L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));
                            %k=1;
                            %Ll1k = double(L(l1,k)); 
                            %myfun_po_m2 = @(m2) fun_po_m2(T,m2,m1,lambda_s,k,gamma_s,f1,f2,f0,phi1,phi2,sigma_1,sigma_2,w,l1,B_0,K_0,M,Ll1k);
                            %m2 = fsolve(myfun_po_m2,0);
                        end
                        k=1;
                        t1 = double(T/k + 1/lambda_s*log((m2-m1)/((exp(lambda_s*T/k)-1)*L(l1,k)+m2-m1*exp(lambda_s*T/k))));
                        fun = symfun(exp(lambda_s*(t1-t))*f(t),t);
                        fun = subs(fun,lambda_s,lambda_s);
                        fun = subs(fun,t1,t1);
                        intfun = double(int(fun, 0, t1));
                        l2 = double((l1+gamma_s*m1/lambda_s)*exp(lambda_s*t1)-gamma_s*m1/lambda_s+gamma_s*K_0(s_position,1)*intfun);
                    else
                        l2 = input('Enter l2 ...');
                        m1 = input('Enter m1 ...');
                        vvod_m2_T = input('Enter 1 if you want to enter "m2" or 2 if you want to enter "T" ...');
                        if vvod_m2_T==1
                            m2 = input('Enter m2 ...');
                            T = input('Enter T(period) ...');
                            H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));       
                            syms H(t)
                            syms L(l1in,k);
                            H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                            L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));
                        else
                            T = input('Enter T(period) ...');
                            m2 = input('Enter m2 ...');
                            %M = (f1*sin(phi1+sigma_1))/sqrt(lambda_s^2+w^2) + (f2*sin(phi2+sigma_2))/sqrt(lambda_s^2+4*w^2);
                            H = double(gamma_s*K_0(s_position,1)*((f1*sin(phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2)));       
                            syms H(t)
                            syms L(l1in,k);
                            H(t) = gamma_s*K_0(s_position,1)*((f1*sin(t + phi1 + sigma_1))/sqrt(lambda_s^2 + w^2)+(f2*sin(2*t + phi2 + sigma_2))/sqrt(lambda_s^2 + 4*w^2));
                            L(l1in,k) = -(lambda_s*l1in/gamma_s + K_0(s_position,1)*f0)+(lambda_s*(H(w*T/k) - H(0)*exp(lambda_s*T/k)))/(gamma_s*(exp(lambda_s*T/k)-1));
                            %k=1;
                            %Ll1k = double(L(l1,k)); 
                            %myfun_po_m2 = @(m2) fun_po_m2(T,m2,m1,lambda_s,k,gamma_s,f1,f2,f0,phi1,phi2,sigma_1,sigma_2,w,l1,B_0,K_0,M,Ll1k);
                            %m2 = fsolve(myfun_po_m2,0);
                        end
                        k=1;
                        t1 = double(T/k + 1/lambda_s*log((m2-m1)/((exp(lambda_s*T/k)-1)*L(l1,k)+m2-m1*exp(lambda_s*T/k))));
                        fun = symfun(exp(lambda_s*(T/k-t))*f(t),t);
                        fun = subs(fun,lambda_s,lambda_s);
                        fun = subs(fun,T,T);
                        fun = subs(fun,k,k);
                        intfun = double(int(fun, t1, T/k));
                        l1 = double((l2+gamma_s*m2/lambda_s)*exp(lambda_s*(T/k-t1))-gamma_s*m2/lambda_s+gamma_s*K_0(s_position,1)*intfun);
                    end
                end
                %%
                %Проверка условия теоремы 1
                if (m2-m1*exp(lambda_s*T/k)+(exp(lambda_s*T/k)-1)*L(l1,k)>0)&&(m1<L(l1,k))&&(L(l1,k)<m2)&&(t1==double(T/k+(1/lambda_s)*log((m2-m1)/((exp(lambda_s*T/k)-1)*L(l1,k)+m2-m1*exp(lambda_s*T/k)))))
                    disp('Система трансцендентных уравнений имеет решение (t1,T/k), где t1 принадлежит от 0 до T/k');
                else
                    disp('Система трансцендентных уравнений НЕ имеет решение (t1,T/k), где t1 принадлежит от 0 до T/k');    
                end
                %%
                %t1 = double(T/k + 1/S(i,1)*log((m2-m1)/((exp(S(i,1)*T/k)-1)*L(l1,k)+m2-m1*exp(S(i,1)*T/k))));
                %Si1=double(S(i,1));
                %fun = @(t,t1,Si1,f0,f1,f2,w,phi1,phi2)exp(Si1*(t1-t)).*(f0+f1*sin(w*t+phi1)+f2*sin(2*w*t+phi2));
                %f_s_podst = @(t) double(subs(f,t,t));
                %fun = @(t,t1,Si1) exp(Si1*(t1-t))*f(t);
                %fun = symfun(exp(Si1*(t1-t))*f(t),t);
                %fun = subs(fun,Si1,Si1);
                %fun = subs(fun,t1,t1);
                %intfun = double(int(fun, 0, t1));
                %intfun = integral(@(t)fun(t,t1,Si1,f0,f1,f2,w,phi1,phi2),0,t1);
                %intfun = integral(@(t)fun(t,t1,Si1),0,t1);
                %l2 = double((l1+gamma_s*m1/S(i,1))*exp(S(i,1)*t1)-gamma_s*m1/S(i,1)+gamma_s*K_0(i,1)*intfun);
                X1 = zeros(raz, 1); 
                X2 = zeros(raz, 1); 
                %fun = @(t,t1,Sj1,f0,f1,f2,w,phi1,phi2)exp(-Sj1*t).*(f0+f1*sin(w*t+phi1)+f2*sin(2*w*t+phi2));
                %fun = symfun(exp(-1*Sj1*t)*f(t),t);
                %intfun = double(int(fun, 0, t1));
                for j=1:1:raz
                    if j~=s_position 
                        Sj1=double(S(j,1));
                        fun = symfun(exp(-1*Sj1*t)*f(t),t);
                        fun_s_pod = subs(fun,Sj1,Sj1);
                        intfun1 = double(int(fun_s_pod, 0, t1));
                        intfun2 = double(int(fun_s_pod, t1, T/k));
                        %intfun1 = double(integral(@(t)fun(t,t1,Sj1,f0,f1,f2,w,phi1,phi2),0,t1));
                        %intfun2 = double(integral(@(t)fun(t,t1,Sj1,f0,f1,f2,w,phi1,phi2),t1,T/k));
                        As = [double(exp(S(j,1)*t1)),-1;
                              -1,double(exp(S(j,1)*((T/k)-t1)))];
                        Bs = [double((m1/S(j,1))*(1-exp(S(j,1)*t1)) - K_0(j,1)*exp(S(j,1)*t1) * intfun1);
                              double( (m2/S(j,1))*(1-exp(S(j,1)*(T/k - t1))) - K_0(j,1)*exp(S(j,1)*T/k) * intfun2)];
                        X = linsolve(As,Bs);
                        X1(j,1) = double(X(1,1));
                        X2(j,1) = double(X(2,1));
                    else
                        X1(j,1) = double(l1/gamma_s);
                        X2(j,1) = double(l2/gamma_s);
                    end
                end
                Y1 = Smatrix*X1;
                Y2 = Smatrix*X2;
                %disp(Y1);
                %disp(Y2);
                %%
                x_like_t1 = linspace(0,T/k,100);
                y_left_part_eq = zeros(1,100);
                y_right_part_eq = zeros(1,100);
                for i=1:100
                    y_left_part_eq(i) = l2;
                    fun = symfun(exp(lambda_s*(t1-t))*f(t),t);
                    fun = subs(fun,lambda_s,lambda_s);
                    fun = subs(fun,t1,x_like_t1(i));
                    intfun = double(int(fun, 0, x_like_t1(i)));
                    y_right_part_eq(i) = double((l1+gamma_s*m1/lambda_s)*exp(lambda_s*(x_like_t1(i)))-gamma_s*m1/lambda_s+gamma_s*K_0(s_position,1)*intfun);
                end
                figure;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                x_toch_peres = x_like_t1;
                y1_toch_peres = y_left_part_eq;
                y2_toch_peres = y_right_part_eq;
                % находим индексы элементов, где разность меняет знак
                d_toch_peres = y2_toch_peres-y1_toch_peres;
                s_toch_peres = abs(diff(sign(d_toch_peres))); 
                id_toch_peres = find( s_toch_peres>0 );

                % пустые массивы (заготовки)
                t_toch_peres = zeros(size(id_toch_peres));
                f_toch_peres = zeros(size(id_toch_peres));
                for k_toch_peres = 1:length(id_toch_peres)
                  i_toch_peres = id_toch_peres(k_toch_peres); % индексы левых точек
                  % находим коэф-ты прямой:
                  a1_toch_peres = (y1_toch_peres(i_toch_peres+1)-y1_toch_peres(i_toch_peres))/(x_toch_peres(i_toch_peres+1)-x_toch_peres(i_toch_peres));
                  a2_toch_peres = (y2_toch_peres(i_toch_peres+1)-y2_toch_peres(i_toch_peres))/(x_toch_peres(i_toch_peres+1)-x_toch_peres(i_toch_peres));
                  b1_toch_peres = y1_toch_peres(i_toch_peres)-a1_toch_peres*x_toch_peres(i_toch_peres);
                  b2_toch_peres = y2_toch_peres(i_toch_peres)-a2_toch_peres*x_toch_peres(i_toch_peres);
                  % имеем два уравнения:
                  % f = a1*t + b1;
                  % f = a2*t + b2;
                  % получаем систему:
                  % -a1*t + f = b1;
                  % -a2*t + f = b2;
                  % матрица коэф-тов:
                  A_toch_peres = [-a1_toch_peres, 1; 
                       -a2_toch_peres, 1];
                  % столбец левой части: 
                  B_toch_peres = [b1_toch_peres; b2_toch_peres];
                  u_toch_peres = A_toch_peres\B_toch_peres; % решаем систему
                  % сохраняем результат:
                  t_toch_peres(k_toch_peres) = u_toch_peres(1);
                  f_toch_peres(k_toch_peres) = u_toch_peres(2);
                end

                % plot them all
                plot(x_toch_peres, y1_toch_peres,'-b', x_toch_peres,y2_toch_peres,'-r', t_toch_peres,f_toch_peres,'ok','linew',2)
                grid on
                %%
                disp('Точки пересечения:');
                for i=1:size(t_toch_peres,2)
                    str_tochki_peres=sprintf('[%.5f, %.5f]',t_toch_peres(i),f_toch_peres(i));
                    disp(str_tochki_peres);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%
                postr_gr = input('Enter 1 if you want to enter the number of turns or 2 if you want to enter a time interval ...');
                if postr_gr == 1
                    nu = input('Enter "nu" ...');
                    t_b = 0;
                    t_e = T*nu/k;
                else
                    t_b = input('Enter left bound of the time interval ...');
                    t_e = input('Enter right bound of the time interval ...');
                end

                S=double(S);
                Q=@(t,m,S,B_0,K_0,i,f0,f1,f2,w,phi1,phi2) (B_0(i,1)*m +K_0(i,1)*f0)/S(i,1) + ((K_0(i,1)*f1)/(S(i,1)^2 +w^2))*(S(i,1)*sin(w*t+phi1)+w*cos(w*t+phi1))+((K_0(i,1)*f2)/(S(i,1)^2 +4*(w^2)))*(S(i,1)*sin(2*w*t+phi2)+2*w*cos(2*w*t+phi2));
                X_0=[double((exp(S(1,1)*(T/k-t1))*(Q(t1,m2,S,B_0,K_0,1,f0,f1,f2,w,phi1,phi2)-Q(t1,m1,S,B_0,K_0,1,f0,f1,f2,w,phi1,phi2))+ exp(S(1,1)*(T/k))*Q(0,m1,S,B_0,K_0,1,f0,f1,f2,w,phi1,phi2)-Q(T,m2,S,B_0,K_0,1,f0,f1,f2,w,phi1,phi2))*(1-exp(S(1,1)*T/k))^(-1));
                     double((exp(S(2,1)*(T/k-t1))*(Q(t1,m2,S,B_0,K_0,2,f0,f1,f2,w,phi1,phi2)-Q(t1,m1,S,B_0,K_0,2,f0,f1,f2,w,phi1,phi2))+ exp(S(2,1)*(T/k))*Q(0,m1,S,B_0,K_0,2,f0,f1,f2,w,phi1,phi2)-Q(T,m2,S,B_0,K_0,2,f0,f1,f2,w,phi1,phi2))*(1-exp(S(2,1)*T/k))^(-1));
                     double((exp(S(3,1)*(T/k-t1))*(Q(t1,m2,S,B_0,K_0,3,f0,f1,f2,w,phi1,phi2)-Q(t1,m1,S,B_0,K_0,3,f0,f1,f2,w,phi1,phi2))+ exp(S(3,1)*(T/k))*Q(0,m1,S,B_0,K_0,3,f0,f1,f2,w,phi1,phi2)-Q(T,m2,S,B_0,K_0,3,f0,f1,f2,w,phi1,phi2))*(1-exp(S(3,1)*T/k))^(-1))];
                x_1=@(t,m1,S,B_0,K_0,i,f0,f1,f2,w,phi1,phi2,nu,k,T,X_0) exp(S(i,1)*(t-(nu-1)*T/k))*(X_0(i,1)+ Q((nu-1)*T/k,m1,S,B_0,K_0,i,f0,f1,f2,w,phi1,phi2))-Q(t,m1,S,B_0,K_0,i,f0,f1,f2,w,phi1,phi2);
                x_2=@(t,m2,S,B_0,K_0,i,f0,f1,f2,w,phi1,phi2,nu,k,T,t1,X2) exp(S(i,1)*(t-(t1+(nu-1)*T/k))) * (X2(i,1)+ Q((t1+(nu-1)*T/k),m2,S,B_0,K_0,i,f0,f1,f2,w,phi1,phi2))-Q(t,m2,S,B_0,K_0,i,f0,f1,f2,w,phi1,phi2);
                t_per=linspace(t_b,t_e,2000); 
                x=double.empty(0,2000);
                y=double.empty(0,2000);
                z=double.empty(0,2000);

                for i=1:2000
                    nu = fix(abs(double(t_per(i)))/double(T/k));
                    nu=nu+1;
                    if (abs(mod(t_per(i),T/k))<=t1) && ( abs(mod(t_per(i),T/k))>=0)    
                        x(i)= x_1(t_per(i),m1,S,B_0,K_0,1,f0,f1,f2,w,phi1,phi2,nu,1,T,X_0);
                        y(i)= x_1(t_per(i),m1,S,B_0,K_0,2,f0,f1,f2,w,phi1,phi2,nu,1,T,X_0);
                        z(i)= x_1(t_per(i),m1,S,B_0,K_0,3,f0,f1,f2,w,phi1,phi2,nu,1,T,X_0);
                    else
                        x(i)= x_2(t_per(i),m2,S,B_0,K_0,1,f0,f1,f2,w,phi1,phi2,nu,1,T,t1,X2);
                        y(i)= x_2(t_per(i),m2,S,B_0,K_0,2,f0,f1,f2,w,phi1,phi2,nu,1,T,t1,X2);
                        z(i)= x_2(t_per(i),m2,S,B_0,K_0,3,f0,f1,f2,w,phi1,phi2,nu,1,T,t1,X2);
                    end
                end 
                % t=linspace(t1,T,1000); 
                % for i=1001:2000
                %     x(i)= x_2(t(i-1000),m2,S,B_0,K_0,1,f0,f1,f2,w,phi1,phi2,1,1,T,t1,X2);
                %     y(i)= x_2(t(i-1000),m2,S,B_0,K_0,2,f0,f1,f2,w,phi1,phi2,1,1,T,t1,X2);
                %     z(i)= x_2(t(i-1000),m2,S,B_0,K_0,3,f0,f1,f2,w,phi1,phi2,1,1,T,t1,X2);
                % end 
                figure;
                plot3(x,y,z);
                grid on;
                hold on
                plot3(X1(1,1),X1(2,1),X1(3,1),'r.','MarkerSize', 20);
                plot3(X2(1,1),X2(2,1),X2(3,1),'r.','MarkerSize', 20);



                [y_g1,z_g1]=meshgrid(linspace(min(y),max(y),10),linspace(min(z),max(z),10));
                x_g1 = zeros(1,10);
                for i=1:10
                    x_g1(i) = (l1-GAMMA(2,1)*y_g1(i)-GAMMA(3,1)*z_g1(i))/GAMMA(1,1);
                end

                x_g2 = zeros(1,10);
                for i=1:10
                    x_g2(i) = (l2-GAMMA(2,1)*y_g1(i)-GAMMA(3,1)*z_g1(i))/GAMMA(1,1);
                end
                surf(x_g1,y_g1,z_g1,'FaceAlpha',0.2,'EdgeColor','none','FaceColor','#EDB120')
                surf(x_g2,y_g1,z_g1,'FaceAlpha',0.2,'EdgeColor','none','FaceColor','#7E2F8E')
                hold off

                %%
                %f = @(tau,f0,f1,f2,w,phi1,phi2) f0 + f1*sin(w*tau+phi1)+f2*sin(2*w*tau + phi2);
                % func=@(t,tau,t1,f0,f1,f2,w,phi1,phi2,A,B,K,m1) exp(-A*(tau-t))*(B*m1+K*(f0 + f1*sin(w*tau+phi1)+f2*sin(2*w*tau + phi2)));
                % intfun=@(t,tau,t1,f0,f1,f2,w,phi1,phi2,A,B,K,m1) integral(@(tau)func(t,tau,t1,f0,f1,f2,w,phi1,phi2,A,B,K,m1),t1,t,'ArrayValued',true);
                % Yt=@(t,t1,tau,f0,f1,f2,w,phi1,phi2,A,B,K,m1,Y1)exp(A*(t-t1))*Y1 + intfun(t,tau,t1,f0,f1,f2,w,phi1,phi2,A,B,K,m1);
                x_init=double.empty(0,2000);
                y_init=double.empty(0,2000);
                z_init=double.empty(0,2000);
                for i=1:2000
                    X_other=[x(i);y(i);z(i)];
                    Y_other=Smatrix*X_other;
                    x_init(i)=Y_other(1,1);
                    y_init(i)=Y_other(2,1);
                    z_init(i)=Y_other(3,1);
                end 
                figure;
                plot3(x_init,y_init,z_init);
                grid on;
                hold on;
                plot3(Y1(1,1),Y1(2,1),Y1(3,1),'r.','MarkerSize', 10);
                plot3(Y2(1,1),Y2(2,1),Y2(3,1),'r.','MarkerSize', 10);

                [y_g1,z_g1]=meshgrid(linspace(min(y_init),max(y_init),50),linspace(min(z_init),max(z_init),50));
                x_g1 = zeros(1,50);
                for i=1:50
                    x_g1(i) = (l1-GAMMA(2,1)*y_g1(i)-GAMMA(3,1)*z_g1(i))/GAMMA(1,1);
                end

                x_g2 = zeros(1,50);
                for i=1:50
                    x_g2(i) = (l2-GAMMA(2,1)*y_g1(i)-GAMMA(3,1)*z_g1(i))/GAMMA(1,1);
                end
                surf(x_g1,y_g1,z_g1,'FaceAlpha',0.2,'EdgeColor','none','FaceColor','#EDB120')
                surf(x_g2,y_g1,z_g1,'FaceAlpha',0.2,'EdgeColor','none','FaceColor','#7E2F8E')
                hold off

                %%
                disp('Система с матрицей A:')
                disp(A);
                disp('векторами B:')
                disp(B);
                disp('K:')
                disp(K);
                disp('C:')
                disp(C);
                fprintf('f(t)=');
                disp(f);
                fprintf('l_1 =');
                disp(l1);
                fprintf('m_1 =');
                disp(m1);
                fprintf('l_2 =');
                disp(l2);
                fprintf('m_2 =');
                disp(m2);
                fprintf('имеет решение с параметрами t1:');          
                disp(t1);
                fprintf('Tf ='); 
                disp(T/k);
                disp('Y1:');
                disp(Y1);
                disp('Y2:');
                disp(Y2);
                fprintf('вычисленные с точностью:'); 
                disp(accuracy);
            else
                disp('Корни характеристического уровнения не удовлетворяют условиям');
            end
        else
            disp('Не выполняется условие обратимости преобразований');
        end
   else
       disp('Не выполняется условие полной управляемости')
   end
else
    disp('Другой вид функции');
end
