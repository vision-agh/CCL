% clc;
close all;
clearvars;

%% Generacja obrazka testowego
% Problem: Czy można zwiększyć latencję nadawania, ale zwiększyć liczbę przekododwań wyników?

TEST_IMG_SIZE = 200;

% -- OBRAZ TESTOWY RECZNY
% testBin =[
%     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
%     1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0;
%     1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0;
%     1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0;
%     1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0;
%     1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0;
%     1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0];

% testBin =[
%     0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0;
%     0 0 0 0 1 0 0 0 1 0 1 0 1 0 1 0 1;
%     1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1;
%     1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1;
%     1 1 1 1 1 0 1 0 1 0 1 0 1 0 1 0 1;
%     1 0 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0;
%     1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% inputBin = zeros(100);
% inputBin(2:8, 5:size(testBin,2)+4) = testBin;

%
% Obrazki losowe
% inputBin  = rand(TEST_IMG_SIZE) > 0.75;
% init = zeros(TEST_IMG_SIZE);
load('Input03.mat');
% init(1:100,1:100) = inputBin;
% inputBin = init;
% load('Input00.mat');

%inputBin = imread('test_26.png') > 0;
%inputBin = imread('test_1.png') > 0;


[V, H] = size(inputBin);

inputBin(1, :) = 0;
inputBin(V, :) = 0;
inputBin(:, 1) = 0;
inputBin(:, H) = 0;

% figure;
% imshow(inputBin,[]);

% Postac kontekstu
% A  B  C  D  E  F
% G  L1 L2 L3 L4

%% PARAMETRY
N = 4;          % Liczba pikseli w wektorze
MAX_MERGE = 2;  % Maksymalna liczba połączeń w jednym wektorze
MAX_LABELS = 1024;%*1024;

%% Algorytm

nextNewLabel = 1;           %  Nowa etykieta, ktora należy nadać [na razie bez odzyskiwania]

context = zeros(2, N+2);    %  Kontekst dla rozwaznanego piksela (wektor + piksel z kolumna z lewej i prawej

L = zeros(size(inputBin));  %  Etykiety
L_real = zeros(size(inputBin)); % Rzeczywiste wyj�cie

% Inicjalizacja tablicy na przekodowania
lut = zeros(1, MAX_LABELS);

% Inicjalizacja tablicy LUT (w sprzecie trzeba to chyba robić na bierząco)
for i=1:MAX_LABELS
    lut(i) = i;
end

% Stos na lancuchy zlaczen
mergeNumber = 0;
merge_stack = zeros(100, 2);
merge_stack_top = 1;
BreakBin = false;
doStack = 0;

% Glowna petla po obrazie
for jj=2:V       % Zaczynamy od 2 linii
    
    last_merge = [0,0];
    %disp('++++++++');
    for ii=1:N:H % Inna konstrukcja niż MK
        %disp('====');
        
        %         if( jj >= 71 && ii >= 21 && ii < 37)
        if( jj >= 71 && ii >= 20)
            BreakBin = true;
        end
        
        % Utworzenie kontekstu z etykiet
        % Piewszy kontekst w linii to specyficzny przypadek, bo nie ma
        % etykiet "z lewej".
        if ii == 1
            context(1:2,2:N+2) = L(jj-1:jj,ii:ii+N);
            context(1:2,1) = 0;  % Brak etykiet -  ustalamy 0
        else
            if ii > H-4 % Ostatni kontekts też jest specyficzny, bo brakuje danych z prawej
                context = L(jj-1:jj,ii-1:ii+N-1);
                context(1:2,N+2) = 0; % Brak etykiet - ustalamy 0
                context(2, 1) = L_local(k);
            else % Typowy kontekst w~środku
                context = L(jj-1:jj,ii-1:ii+N);
                context(2, 1) = L_local(k);
            end
        end
        
        % Wektor danych wejściowych (binarny) -- N - pikseli
        P = inputBin(jj,ii:ii+N-1);
        
        % Przepuszczenie kontekstu przez LUT tj. nadawanie najnowszych
        % etykiet
        % W HW trzeba zrobić tak jak PC -jeśli zmiana nastąpiła w
        % poprzednim takcie, to potrzebne są specjalnie multipleksery
        for i=1:N+2
            % Przekodowanie wyjścia z linii opóźniających
            if ( context(1,i) > 0 )
                context(1,i) = lutD(context(1,i));
                for merge_idx = 1:merge_i-1
                    if context(1, i) == merge(merge_idx, 1)
                        context(1, i) = merge(merge_idx, 2);
                        break;
                    end
                end
            end
        end
        
        % Przekodowanie ostatnio nadanej etykiety [pytanie czy to powinno się tutuaj dziać]
        if (context(2,1) > 0)
            for merge_idx = 1:merge_i-1
                if context(2, 1) == merge(merge_idx, 1)
                    context(2, 1) = merge(merge_idx, 2);
                    break;
                end
            end
        end
        
        
        % Przypisanie pikseli
        % Może da się taką analizę w ramach jednego bloku Verilog - do przemyslenia/sprawdzenia
        L_local = zeros(1,N);  % Etykiety nadane w ramach kontekstu
        
        
        merge = ones(N,2)*-1;   % Pamieć podręczna na merge
        stack = zeros(N,1);
        merge_i = 1;            % Wskazanik na strukture merge
        if (merge_i > 3)
            breakBin = 1;
        end
        
        % Petla po pikselach w wektorze [w Verilog zrownoleglona]
        for k=1:N
            if ( P(k) == 1 ) % Jesli analizowany piksel nalezy do obiektu
                % Rozpatrujemy 4 podkonteksty 3+1
                
                if ( k == 1 ) % Pierwszy piksel z kontektstu [inny s�siad z lewej]
                    G = context(2,1);
                else
                    G = L_local(k-1);
                end
                
                topRow = context(1, k:k+2);
                switchBin = [G topRow] ~= 0;
                doMerge = 0;
%                 doStack = 0;
                
                switch dec2bin(switchBin)'
                    case '1111'
                        L_local(k) = G;
                    case '1101'
                        if G ~= context(1,k+2)
                            [L_min, L_max] = MinMax(G, context(1,k+2));
                            L_local(k) = L_min;
                            doMerge = 1;
                            doStack = G > context(1,k+2);
                        else
                            L_local(k) = G;
                            doMerge = 0;
%                             doStack = 0;
                        end
                    case '1011'
                        L_local(k) = G;
                    case '1110'
                        L_local(k) = G;
                    case '1100'
                        L_local(k) = G;
                    case '1010'
                        L_local(k) = G;
                    case '1001'
                        if G ~= context(1,k+2)
                            [L_min, L_max] = MinMax(G, context(1,k+2));
                            L_local(k) = L_min;
                            doMerge = 1;
                            doStack = G > context(1,k+2);
                        else
                            L_local(k) = G;
                            doMerge = 0;
%                             doStack = 0;
                        end
                    case '1000'
                        L_local(k) = G;
                    case '0111'
                        L_local(k) = context(1,k+1);
                    case '0101'
                        if context(1,k) ~= context(1,k+2)
                            [L_min, L_max] = MinMax(context(1,k), context(1,k+2));
                            L_local(k) = L_min;
                            doMerge = 1;
                            doStack = context(1,k) > context(1,k+2);
                        else
                            L_local(k) = context(1,k);
                            doMerge = 0;
%                             doStack = 0;
                        end
                    case '0011'
                        L_local(k) = context(1,k+1);
                    case '0110'
                        L_local(k) = context(1,k+1);
                    case '0100'
                        L_local(k) = context(1,k);
                    case '0010'
                        L_local(k) = context(1,k+1);
                    case '0001'
                        L_local(k) = context(1,k+2);
                    case '0000'
                        L_local(k) = nextNewLabel;
                        nextNewLabel = nextNewLabel +1;
                end
                
                if doMerge
                    merge(merge_i,:) = [L_max, L_min];  % Zapisujemy etykiety mniejsz� i wi�ksz�
                    merge_i = merge_i + 1;                   % Inkrementujemy wskaznik merge
                end
                
%                 if doStack
%                     merge_stack(merge_stack_top,:) = [L_max L_min];
%                     merge_stack_top = merge_stack_top +1;
%                     chainFlag = true;
%                 else
%                     chainFlag = false;
%                 end
                
            end
        end
        
        if(BreakBin)
            BreakBin = 1;
        end
        
        if merge_i > 2
            if merge(1, 1) == merge(2, 1)
                if merge(1, 2) > merge(2, 2)
                    merge(1, 1) = merge(1, 2);
                    merge(1, 2) = merge(2, 2);
                else
                    merge(2, 1) = merge(2, 2);
                    merge(2, 2) = merge(1, 2);
                end
                stack = [1; 1];%[1; 1];
            elseif merge(1, 1) == merge(2, 2)
                merge(2, 2) = merge(1, 2);
                stack = [0; 0];
            elseif merge(1, 2) == merge(2, 1)
                merge(1, 2) = merge(2, 2);
                stack = [1; 1];
            else
                stack = [1; 0];
            end
        else
            stack(1) = doStack;
        end
        
        if(BreakBin)
            BreakBin = 1;
        end
        
        % Obsluga merge etykiet [w zależności od liczby tych etykiet może nastąpić wstrzymanie potoku]
        L_localLUT = L_local;
        if (BreakBin)
            BreakBin2 = true;
        end
        lutD = lut;
        merge_i_t = 1;  % licznik tymczasowym
        while merge_i_t < merge_i
            mergeNumber = mergeNumber + 1;
            
            if merge(merge_i_t,1) ~= merge(merge_i_t,2)

                lut(merge(merge_i_t,1)) = merge(merge_i_t,2);   % Zapis połączenia do lut
                mergeP = merge(merge_i_t,1);
                mergeS = merge(merge_i_t,2);

                if stack(merge_i_t)
                    merge_stack(merge_stack_top,:) = merge(merge_i_t, :);
                    merge_stack_top = merge_stack_top + 1;
                end

                %end
                % debug info
                %             if merge_i_t == 1
                %                 X = sprintf('MERGE %d %d',merge(merge_i_t,1),merge(merge_i_t,2));
                %             else
                %                 X = sprintf('DMERGE %d %d',merge(merge_i_t,1),merge(merge_i_t,2));
                %             end
                %             disp(X);
                
                for k = 1:N
                    if ( L_localLUT(k) == mergeP)
                        % Pr�ba. Bez tego powinno dzia�a�.
%                         L(jj,ii+k-1) = lut(L_local(k));
                        L_localLUT(k) = mergeS;
                    end
                end
            end
            
            merge_i_t = merge_i_t + 1;
            
        end
        
        % Zapis etykiet wynikowych
        % Powinniśmy przekazać dalej/zapisać już nowe etykiety - po merge
        for k = 1:N
            if ( L_localLUT(k) > 0 )
%                 L(jj,ii+k-1) = lut(L_local(k));
                L(jj,ii+k-1) = L_localLUT(k);
                L_real(jj,ii+k-1) = L_localLUT(k);
            end
        end
    end
    
    % Operacje po koncu linii [w czasie wygaszania]
    if (merge_stack_top > 1 )
        % E Sklejona | E nadana
        %TODO Ten dodatkowy warunek opisany u Ciaracha [ale to jest optymalizcja]
        while ( merge_stack_top >  1)
            merge_stack_top = merge_stack_top - 1;
            lut(merge_stack(merge_stack_top,1)) = lut(merge_stack(merge_stack_top,2));
            %             X = sprintf('MERGE CHAIN %d %d',merge_stack(merge_stack_top,1),lut(merge_stack(merge_stack_top,2)));
            %             disp(X);
        end
    end
    lutD = lut;
end

%return
L_test = L_real;

% Korekta tablicy sklejen ?
%TODO CZY TO WYSTARCZY (może się trafić case, że nie ?)
for i=1:numel(lut)
    lut(i) = lut(lut(i));
end

% lut = zeros(1,1024);
% for i=1:1024
%     lut(i)=i;
% end



%% Końcowe przekodowanie testowe
for jj=2:V       % Zaczynamy od 2 linii
    for ii=1:1:H
        if ( L(jj,ii) > 0 )
            L(jj,ii) = lut(L(jj,ii));
        end
    end
end

Ln_FPGA = numel(unique(L(:)))-1;
[L_M,Ln_M] = bwlabel(inputBin);

if (Ln_FPGA == Ln_M)
%     disp('FPGA == MATLAB !!!');

else
    disp('FPGA ERROR !!!');
end


%% Porowananie działania z Matlabem
for l=1:Ln_M
    % Pobierz etykiete
    B = (L_M == l);
    % Wytnij z obrazu 4K
    C = L.*B;
    % Sprawdz, czy obiekt ma tylko jedną etykietę
    U = unique(C(:));
    if ( numel(U) > 2 )
        disp('FPGA ERROR !!!')
        U
        %         imwrite(uint8((C>0)*255),'test_2.png');
        %         imwrite(uint8(inputBin*255),'test_1.png');
        
        %         figure(2);
        %         imshow(C,[]);
        
    end
    
    %pause;
end



