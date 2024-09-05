%% Plot results
if strcmpi(paramName2,'M1')
    paramName2 = 'spatial oversampling rate';
    paramSet2 = paramSet2 / N;
end

%% create string for legend
legend_str = {};
for iAlgo = 1:nAlgo
    for iParam1 = 1:nParam1
        legend_str = [legend_str; {sprintf('%s, %s=%d',func2str(algoSet{iAlgo}),paramName1,paramSet1(iParam1))}];
    end
end

%% Complexity
figure;
legend_str = {};
for iAlgo = 2:nAlgo
    semilogx(paramSet2,sum(avgResults.complexity.Niter(:,:,iAlgo,:),4).');
    hold on;
    for iParam1 = 1:nParam1
        legend_str = [legend_str; {sprintf('%s, %s=%d',func2str(algoSet{iAlgo}),paramName1,paramSet1(iParam1))}];
    end
end
hold off;
grid on;
xlabel(paramName2);
ylabel('number of iterations');
legend(legend_str,'Interpreter','none');
matlab2tikz(sprintf('Niter_vary%s.tikz',paramName2));


figure;
legend_str = {};
for iAlgo = [1,3]
    semilogx(paramSet2,sum(avgResults.complexity.runtime(:,:,iAlgo,:),4).');
    hold on;
    for iParam1 = 1:nParam1
        legend_str = [legend_str; {sprintf('%s, %s=%d',func2str(algoSet{iAlgo}),paramName1,paramSet1(iParam1))}];
    end
end
hold off;
grid on;
xlabel(paramName2);
ylabel('CPU time (seconds)');
legend(legend_str,'Interpreter','none');
matlab2tikz(sprintf('cputime_vary%s.tikz',paramName2));


%% Fmeasure
figure;
legend_str = {};
for iAlgo = 1:nAlgo
    semilogx(paramSet2,sum(avgResults.pattern.FmeasureZ(:,:,iAlgo,:),4).');
    hold on;
    for iParam1 = 1:nParam1
        legend_str = [legend_str; {sprintf('%s, %s=%d',func2str(algoSet{iAlgo}),paramName1,paramSet1(iParam1))}];
    end
end
semilogx(paramSet2,avgResults.pattern.FmeasureZ_ref','-.');
hold off;
grid on;
xlabel(paramName2);
ylabel('F-measure(Z)');
legend(legend_str,'Interpreter','none');
matlab2tikz(sprintf('FmeasureZ_vary%s.tikz',paramName2));

% %% Error
% figure;
% plot(avgResults.error.X(:,:,1)');
% hold on;
% plot(avgResults.error.DZ(:,:)','--');
% hold on;
% plot(avgResults.error.X(:,:,2)','-.');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorX');
% % matlab2tikz('errorX.tikz');
% 
% figure;
% plot(avgResults.error.D(:,:,1)');
% hold on;
% plot(avgResults.error.D(:,:,2)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorD');
% % matlab2tikz('errorD.tikz');
% 
% figure;
% plot(avgResults.error.Z(:,:,1)');
% hold on;
% plot(avgResults.error.Z(:,:,2)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorZ');
% % matlab2tikz('errorZ.tikz');
% 
%% Debiased Error
% figure;
% semilogx(paramSet2,20*log10(avgResults.debiasedError.X(:,:,1)).');
% hold on;
% semilogx(paramSet2,20*log10(avgResults.debiasedError.DZ(:,:)).','--');
% hold on;
% semilogx(paramSet2,20*log10(avgResults.debiasedError.X(:,:,2)).','-.');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('MNSE(X)');
% legend(legend_str,'Interpreter','none');
% matlab2tikz(sprintf('errorX_vary%s.tikz',paramName2));

figure;
legend_str = {};
for iAlgo = 2:3
    semilogx(paramSet2,20*log10(avgResults.debiasedError.D(:,:,iAlgo)).');
    hold on;
    for iParam1 = 1:nParam1
        legend_str = [legend_str; {sprintf('%s, %s=%d',func2str(algoSet{iAlgo}),paramName1,paramSet1(iParam1))}];
    end
end
hold off;
grid on;
xlabel(paramName2);
ylabel('MNSE(D)');
legend(legend_str,'Interpreter','none');
matlab2tikz(sprintf('errorD_vary%s.tikz',paramName2));

figure;
legend_str = {};
for iAlgo = [2,3]
    semilogx(paramSet2,20*log10(avgResults.debiasedError.Z(:,:,iAlgo)).');
    hold on;
    for iParam1 = 1:nParam1
        legend_str = [legend_str; {sprintf('%s, %s=%d',func2str(algoSet{iAlgo}),paramName1,paramSet1(iParam1))}];
    end
end
hold off;
grid on;
xlabel(paramName2);
ylabel('MNSE(Z)');
legend(legend_str,'Interpreter','none');
matlab2tikz(sprintf('errorZ_vary%s.tikz',paramName2));


% 
% %% Error before debias vs after debias for formulation w\o X
% figure;
% plot(avgResults.debiasedError.X(:,:,1)');
% hold on;
% plot(avgResults.error.X(:,:,1)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorX');
% title('without X');
% % matlab2tikz('errorX_debias.tikz');
% 
% figure;
% plot(avgResults.debiasedError.D(:,:,1)');
% hold on;
% plot(avgResults.error.D(:,:,1)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorD');
% title('without X');
% % matlab2tikz('errorD_debias.tikz');
% 
% figure;
% plot(avgResults.debiasedError.Z(:,:,1)');
% hold on;
% plot(avgResults.error.Z(:,:,1)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorZ');
% title('without X');
% matlab2tikz('errorZ_debias.tikz');
% 
% %% Error before debias vs after debias for formulation w\ X
% figure;
% plot(avgResults.debiasedError.X(:,:,2)');
% hold on;
% plot(avgResults.error.X(:,:,2)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorX');
% title('with X');
% % matlab2tikz('errorX_debias.tikz');
% 
% figure;
% plot(avgResults.debiasedError.DZ(:,:)');
% hold on;
% plot(avgResults.error.DZ(:,:)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorDZ');
% title('with X');
% % matlab2tikz('errorX_debias.tikz');
% 
% figure;
% plot(avgResults.debiasedError.D(:,:,2)');
% hold on;
% plot(avgResults.error.D(:,:,2)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorD');
% title('with X');
% % matlab2tikz('errorD_debias.tikz');
% 
% figure;
% plot(avgResults.debiasedError.Z(:,:,2)');
% hold on;
% plot(avgResults.error.Z(:,:,2)','--');
% hold off;
% grid on;
% xlabel(paramName2);
% ylabel('errorZ');
% title('with X');
% % matlab2tikz('errorZ_debias.tikz');