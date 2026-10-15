T = readtable('big12.csv', 'VariableNamingRule','preserve');

season    = T.season;
week      = T.week;
gameDate  = datetime(T.date);
home      = string(T.home);
home_pts  = T.home_pts;
away      = string(T.away);
away_pts  = T.away_pts;
completed = T.completed;

teams = unique([home; away], 'sorted');
m = numel(teams);

fprintf('Number of teams: %d\n',m);

H25 = buildH(teams,season,home,away,home_pts,away_pts,...
             completed,2025);
H26 = buildH(teams,season,home,away,home_pts,away_pts,...
             completed,2026);

PiFlat = 0.5*ones(m);
PiFlat(1:m+1:end) = 0;

PiLast = H25;
missing = isnan(PiLast);
PiLast(missing) = 0.5;
PiLast(1:m+1:end) = 0;

strength = ones(m,1);
PiPre = zeros(m);

for i = 1:m
    for j = 1:m
        if i ~= j
            PiPre(i,j) = strength(i)/(strength(i)+strength(j));
        end
    end
end

priors = {PiFlat,PiLast,PiPre};
names = {'Flat','2025 head-to-head','Preseason'};

w = 1;
results = struct();

figure('Color','w','Position',[100 100 1100 850]);

for p = 1:3
    Pi = priors{p};

    A = H26 + w*Pi;
    A(isnan(A)) = 0;

    [V,D] = eig(A);
    ev = diag(D);
    [~,order] = sort(abs(ev),'descend');
    ev = ev(order);
    V = V(:,order);

    lambdaStar = ev(1);
    lambda2 = ev(2);

    r = real(V(:,1));
    if sum(r) < 0
        r = -r;
    end
    r = abs(r);
    r = r/sum(r);

    results(p).A = A;
    results(p).r = r;
    results(p).lambdaStar = lambdaStar;
    results(p).lambda2 = lambda2;
    results(p).ratio = abs(lambda2/lambdaStar);

    [rSorted,idx] = sort(r,'descend');

    fprintf('\n===== %s prior =====\n',names{p});
    fprintf('lambda* = %.8f\n',real(lambdaStar));
    fprintf('|lambda2/lambda*| = %.8f\n',...
        abs(lambda2/lambdaStar));

    disp(table(teams(idx),rSorted,...
        'VariableNames',{'Team','RankingScore'}));

    subplot(3,1,p);
    bar(rSorted);
    xticks(1:m);
    xticklabels(teams(idx));
    xtickangle(45);
    ylabel('Ranking score r_i');
    title(names{p});
    grid on;
end

sgtitle('2026 Big 12 rankings, w = 1');
exportgraphics(gcf,'big12_rankings.png','Resolution',200);

chosenPrior = 1;
r = results(chosenPrior).r;

targetDates = datetime(2026,10,[10 11]);
gameMask = season == 2026 & ...
           dateshift(gameDate,'start','day') >= targetDates(1) & ...
           dateshift(gameDate,'start','day') <= targetDates(2);

gameRows = find(gameMask);

fprintf('\nPredictions for October 10-11, 2026\n');
fprintf('Model: %s prior\n',names{chosenPrior});

correct = 0;
evaluated = 0;

for k = 1:numel(gameRows)
    g = gameRows(k);

    hi = find(teams == home(g),1);
    ai = find(teams == away(g),1);

    if r(hi) > r(ai)
        pred = home(g);
    elseif r(ai) > r(hi)
        pred = away(g);
    else
        pred = "Tie";
    end

    fprintf('%s vs %s: predicted winner = %s\n',...
        home(g),away(g),pred);

    if completed(g) == 1
        if home_pts(g) > away_pts(g)
            actual = home(g);
        elseif away_pts(g) > home_pts(g)
            actual = away(g);
        else
            actual = "Tie";
        end

        isCorrect = (pred == actual);
        correct = correct + isCorrect;
        evaluated = evaluated + 1;

        fprintf('  Final: %g-%g; actual winner: %s; correct: %d\n',...
            home_pts(g),away_pts(g),actual,isCorrect);
    else
        fprintf('  Game not completed; not scored yet.\n');
    end
end

if evaluated > 0
    fprintf('\nAccuracy: %d/%d = %.1f%%\n',...
        correct,evaluated,100*correct/evaluated);
else
    fprintf('\nNo completed games found in the date range.\n');
end

function H = buildH(teams,season,home,away,home_pts,away_pts,...
                    completed,targetSeason)

    m = numel(teams);
    H = nan(m,m);
    H(1:m+1:end) = 0;

    sumH = zeros(m,m);
    countH = zeros(m,m);

    rows = find(season == targetSeason & completed == 1);

    for k = 1:numel(rows)
        g = rows(k);
        i = find(teams == home(g),1);
        j = find(teams == away(g),1);

        hp = double(home_pts(g));
        ap = double(away_pts(g));

        if hp + ap <= 0
            continue;
        end

        sumH(i,j) = sumH(i,j) + hp/(hp+ap);
        sumH(j,i) = sumH(j,i) + ap/(hp+ap);

        countH(i,j) = countH(i,j) + 1;
        countH(j,i) = countH(j,i) + 1;
    end

    for i = 1:m
        for j = 1:m
            if i ~= j && countH(i,j) > 0
                H(i,j) = sumH(i,j)/countH(i,j);
            end
        end
    end
end
