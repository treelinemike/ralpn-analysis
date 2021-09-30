% display scope tracker visibility
% Requires GitHub surgnav-tools repo in MATLAB path

% restart
close all; clear; clc;

% options
doSaveVizPlots = 1;
doMakeVideo = 1;
dt = 0.1;            % [sec] resampling time interval
movAvgWindow = 2*60; % [sec] width of moving average filter
true_fields = {'case_id','gender','laterality','approach','track_file','scope_tool_id','t_end_exp','t_end_loc','t_end_isc','t_undock'};
base_path = 'S:/eit/mak/RALPN-D19002-31584';
% base_path = 'C:/Users/f002r5k/Dropbox/projects/surg_nav/nccc_pilot';
param_filenames = { ...
    '/31584-001/tracker/31584-001-tracker-params.txt';
    '/31584-002/tracker/31584-002-tracker-params.txt';
    '/31584-003/tracker/31584-003-tracker-params.txt';
    '/31584-004/tracker/31584-004-tracker-params.txt';
    '/31584-005/tracker/31584-005-tracker-params.txt';
    '/31584-006/tracker/31584-006-tracker-params.txt';
    '/31584-007/tracker/31584-007-tracker-params.txt';
    '/31584-008/tracker/31584-008-tracker-params.txt';
    '/31584-009/tracker/31584-009-tracker-params.txt';
    '/31584-010/tracker/31584-010-tracker-params.txt';
    '/31584-011/tracker/31584-011-tracker-params.txt';
    '/31584-012/tracker/31584-012-tracker-params.txt';
    };

% define phases and set colors
phaseStrings = {'Exposure', 'Tumor Localization', 'Warm Ischemia','Renorrhaphy & Closing'};
colors = getCustomColors;
phaseColors = [ colors(3:4,:); [1 0 0]; colors(8,:)];

% storage variables for
phaseViz = zeros(length(param_filenames),length(phaseStrings));
phaseDur = zeros(length(param_filenames),length(phaseStrings));

% iterate through each case
for caseIdx = 1:length(param_filenames)
    %% load tracking parameters from file
    thisFileFullPath = [base_path param_filenames{caseIdx}];
    disp(['Opening ' thisFileFullPath]);
    fid = fopen(thisFileFullPath);
    paramArray = textscan(fid,'%s %s','Delimiter',{'=',' '},'MultipleDelimsAsOne',1,'CollectOutput',1);
    fclose(fid);
    if(~prod(size(paramArray{1}) == [length(true_fields),2]))
        error('Incorrect parameter array size');
    end
    paramArray = paramArray{1};
    params = [];
    for keyIdx = 1:size(paramArray,1)
        thisVal = str2double(paramArray{keyIdx,2});
        if(isnan(thisVal))
            thisVal =  paramArray{keyIdx,2};
        end
        params.(paramArray{keyIdx,1}) = thisVal;
    end
    
    % validate keys in params struct
    params_fields = fieldnames(params);
    
    for fieldIdx = 1:length(true_fields)
        if( isempty( find(strcmp(true_fields{fieldIdx},params_fields))   ))
            error(['Parameter file does not include "' true_fields{fieldIdx} '" key']);
        end
    end
    
    % add path to tracker file
    [mat,tok] = regexp(param_filenames{caseIdx},'(\S+[\\/])\S+$','match','tokens');
    if( ~isempty(tok{1}) )
        params.track_file = [base_path tok{1}{1} params.track_file];
    end
    
    %% load actual tracker data
    % load tracker data
    thisFileFullPath = params.track_file;
    disp(['Opening ' thisFileFullPath]);
    fid = fopen(thisFileFullPath);
    allData = textscan(fid,'%f %u %c %f %f %f %f %f %f %f %*f %s','Delimiter',',');
    fclose(fid);
    time = allData{1};
    toolID = allData{3};
    Q = [allData{4:7}];
    T = [allData{8:10}];
    trackLabels = allData{11};
    
    %% compute resampled scope visibility
    scopeTime = time(toolID == params.scope_tool_id); % [sec]
    scopeTimeMin = scopeTime/60;  % [min]
    stdTimes = 0:dt:(params.t_undock*60); % [sec]
    disp('Interpolating scope visibility, may take several seconds...');
    tic
    scopeViz = any((scopeTime > (stdTimes - (dt/2))) & (scopeTime < (stdTimes + (dt/2)) )); % this is VERY slow
    toc
    
    %% average visibility for specified phases
    phaseEndpts = [0, params.t_end_exp, params.t_end_loc, params.t_end_isc, params.t_undock]*60; % note: converted to seconds!
    for thisEndpt = 2:length(phaseEndpts)
        phaseStartIdx = find(stdTimes > phaseEndpts(thisEndpt-1),1,'first');
        phaseEndIdx   = find(stdTimes < phaseEndpts(thisEndpt),1,'last');
        thisPhaseViz  = scopeViz(phaseStartIdx:phaseEndIdx);
        phaseViz(caseIdx,thisEndpt-1) = sum(thisPhaseViz)/length(thisPhaseViz);
    end
    phaseDur(caseIdx,:) = diff(phaseEndpts)/60; % [min]
    
    %% plot results
    figure;
    set(gcf,'Position',[0105 0477 5.736000e+02 0285]);
    hold on; grid on;
    ylim([-0.25,1.125]);
    xlim([0,max(scopeTimeMin)]);
    yticklabels({'0%','50%','100%'});
    yticks([0 0.5 1]);
    
    % plot scope visibility along with moving average filtered signal
    % plot([0;scopeTimeMin],[0;gradient(scopeTimeMin)<0.15],'b-');
    % plot(stdTimes/60,scopeViz,'b-');
    % area(stdTimes/60,scopeViz,'FaceColor',[0.6 0.6 0.6],'LineStyle','none');
    
    bar(stdTimes/60,scopeViz,1,'FaceColor',[0.7 0.7 0.7],'LineStyle','none');
    plot(stdTimes/60,movmean(scopeViz,movAvgWindow/dt),'k-','LineWidth',1.6)
    
    % plot phase durations below scope visibility
    legPlots = [];
    legPlots(end+1) = plot([0 params.t_end_exp],-0.125*ones(1,2),'-','LineWidth',4,'Color',phaseColors(1,:));
    legPlots(end+1) = plot([params.t_end_exp params.t_end_loc],-0.125*ones(1,2),'-','LineWidth',4,'Color',phaseColors(2,:));
    legPlots(end+1) = plot([params.t_end_loc params.t_end_isc],-0.125*ones(1,2),'-','LineWidth',4,'Color',phaseColors(3,:));
    legPlots(end+1) = plot([params.t_end_isc params.t_undock],-0.125*ones(1,2),'k-','LineWidth',4,'Color',phaseColors(4,:));
    
    % generate legend labels, replacing first space in phase strings with
    % newline character
    legStrings = {};
    for phaseStringIdx = 1:length(phaseStrings)
        thisPhaseString = phaseStrings{phaseStringIdx};
        firstSpaceIdx = find(phaseStrings{phaseStringIdx} == ' ',1,'first');
        if(~isempty(firstSpaceIdx))
            thisPhaseString(firstSpaceIdx) = char(10);
        end
        legStrings{end+1} = thisPhaseString;
    end
    
    % extract subject number
    [mat,tok] = regexp(params.case_id,'\d+-(\d+)','match','tokens');
    subjNum = str2num(tok{1}{1});
    
    % generate legend and finalize plot
    legh = legend(legPlots,legStrings,'Orientation','horizontal','Location','southoutside');
    xlabel(['\bfElapsed Time [min]']);
    title(['\bf Endoscope Visibility: ' sprintf('%1.0f',movAvgWindow/60) 'min Moving Average']);
    %     title(['\bfEndoscope Tracker Visibility - ' params.case_id],'FontSize',12);
    %      title(['\bf' params.case_id ': ' params.laterality ' ' params.approach],'FontSize',12);
    title(sprintf('Subject %03d: %s %s',subjNum,params.laterality,params.approach));
    ylabel('\bfScope Visibility');
    xlim([0,200]);
    drawnow;
    
    % save plot if desired
    if(doSaveVizPlots)
        thisImgFile = sprintf('frame%03d.png',caseIdx);
        saveas(gcf,thisImgFile);
        system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
    end
    
end


% compile animation
if(doSaveVizPlots && doMakeVideo)
    system(['ffmpeg -y -r 1 -start_number 1 -i frame%003d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 scope_viz.mp4']);
%     system('del frame*.png');
end

% % % phaseVizA = [
% % %     0.4313    0.3216    0.0623    0.0125
% % %     0.1168    0.0002    0.0001    0.0085
% % %     0.4499    0.5865    0.5361    0.2412
% % %     0.5757    0.3447    0.3474    0.1669
% % %     0.3378    0.6265    0.9400    0.7551];
% % % phaseDurA = [
% % %     72.0000   11.0000   15.0000   42.6150
% % %    54.8300    9.0700   24.0500   34.0836
% % %    88.2200   13.9000   18.8600   23.6442
% % %    61.7300   23.8700   10.4500   51.8572
% % %    26.6200    9.4800   16.7500   48.7872];

%% statistics

% duration of exposure phase
x = phaseDur(:,1);
sem_x = std(x)/sqrt(length(x));
ci_low = mean(x) -sem_x*tinv(0.975,length(x)-1);
ci_high = mean(x) + sem_x*tinv(0.975,length(x)-1);
fprintf('Exposure duration: %03.1f min (95%% CI: %03.1f, %03.1f)\n',mean(x),ci_low,ci_high);

% overall duration of robotic part of procedure
x = sum(phaseDur,2);
sem_x = std(x)/sqrt(length(x));
ci_low = mean(x) -sem_x*tinv(0.975,length(x)-1);
ci_high = mean(x) + sem_x*tinv(0.975,length(x)-1);
fprintf('Mean Robot Minutes: %03.1f min (95%% CI: %03.1f, %03.1f)\n',mean(x),ci_low,ci_high);
fprintf('Robot Minutes Median & Range: %03.1f min (%03.1f, %03.1f)\n',median(x),min(x),max(x));

x = phaseViz(:,1);
sem_x = std(x)/sqrt(length(x));
ci_low = mean(x) -sem_x*tinv(0.975,length(x)-1);
ci_high = mean(x) + sem_x*tinv(0.975,length(x)-1);
fprintf('Exposure scope visibility: %03.1f%% (95%% CI: %03.1f, %03.1f)\n',mean(x)*100,ci_low*100,ci_high*100);


%% Composite Bar Charts

% Composite Tracker Visibility
figure;
% set(gcf,'Position',[1.266000e+02 0342 9.214000e+02 3.742000e+02]);
set(gcf,'Position',[1.266000e+02 4.698000e+02 9.214000e+02 2.464000e+02]);
hold on; grid on;
% phaseStringsTemp = phaseStrings;
% for i = 1:length(phaseStringsTemp)
%     phaseStringsTemp{i} = sprintf('%s \n test',phaseStringsTemp{i});
% end
x = categorical(phaseStrings,phaseStrings);
xbar = mean(phaseViz*100); % [percentage]
sem = std(phaseViz*100)/sqrt(size(phaseViz,1)); % [percentage]
bar(x,xbar,'CData',phaseColors,'FaceColor','flat','EdgeColor','k','LineWidth',1.0);
erb = errorbar(x,xbar,sem,sem);
erb.Color = 'k';
erb.LineStyle = 'none';
erb.LineWidth = 1.0;
title(sprintf('\\bfMean Endoscope Tracker Visibility by Phase (n=%1d)',length(param_filenames)),'FontSize',12);
ylabel('\bfMean Visibility');
ytickformat(gca, 'percentage');
ylim([0,100]);
ax = gca;
ax.LineWidth = 1.0;
ax.XTickLabel = cellfun(@(x) ['\bf' x], ax.XTickLabel, 'UniformOutput', false);
for i = 1:length(x)
    text(x(i),xbar(i)-sem(i)-8*max(get(gca,'YLim'))/100,sprintf('\\bf%0.1f \x00B1 %0.1f',xbar(i),sem(i)),'HorizontalAlignment','Center');
end

% Composite Phase Duration
figure;
% set(gcf,'Position',[1.266000e+02 0342 9.214000e+02 3.742000e+02]);
set(gcf,'Position',[1.266000e+02 4.698000e+02 9.214000e+02 2.464000e+02]);
hold on; grid on;
x = categorical(legStrings,legStrings);
% yyaxis right;
% phaseFrac = phaseDur ./ repmat(sum(phaseDur,2),1,4);
% xbar2 = mean(phaseFrac*100); % [percentage]
% bar(x,xbar2,'CData',phaseColors,'FaceColor','flat','EdgeColor','k','LineWidth',1.0);
% sem = std(phaseFrac*100)/sqrt(size(phaseViz,1)); % [percentage]
% erb = errorbar(x,xbar,sem,sem);
% erb.Color = 'k';
% erb.LineStyle = 'none';
% erb.LineWidth = 1.0;
% ylabel('\bfFraction of Total Robot Time');
% ylim([0,100*80/sum(mean(phaseDur))]);
% ax = gca;
% ytickformat(ax, 'percentage');
% ax.YColor = [0 0 0];
%
% yyaxis left;
xbar = mean(phaseDur); % [min]
sem = std(phaseDur)/sqrt(size(phaseViz,1)); % [min]
bar(x,xbar,'CData',phaseColors,'FaceColor','flat','EdgeColor','k','LineWidth',1.0);
erb = errorbar(x,xbar,sem,sem);
erb.Color = 'k';
erb.LineStyle = 'none';
erb.LineWidth = 1.0;

title(sprintf('\\bfMean Phase Duration (n=%1d)',length(param_filenames)),'FontSize',12);
ylabel('\bfMean Duration [min]');
ylim([0 80]);
ax = gca;
ax.YColor = [0 0 0];
ax.LineWidth = 1.0;
ax.XTickLabel = cellfun(@(x) ['\bf' x], ax.XTickLabel, 'UniformOutput', false);
for i = 1:length(x)
    text(x(i),xbar(i)-sem(i)-8*max(get(gca,'YLim'))/100,sprintf('\\bf%0.1f \x00B1 %0.1f',xbar(i),sem(i)),'HorizontalAlignment','Center');
end
% Composite Phase Fraction
phaseFrac = phaseDur ./ repmat(sum(phaseDur,2),1,4);
figure;
hold on; grid on;
x = categorical(phaseStrings,phaseStrings);
xbar = mean(phaseFrac*100); % [percentage]`
sem = std(phaseFrac*100)/sqrt(size(phaseViz,1)); % [percentage]
bar(x,xbar,'CData',phaseColors,'FaceColor','flat','EdgeColor','k','LineWidth',1.0);
erb = errorbar(x,xbar,sem,sem);
erb.Color = 'k';
erb.LineStyle = 'none';
erb.LineWidth = 1.0;
ylabel('\bfFraction of Total Robot Time');
ytickformat(gca, 'percentage');
ax3 = gca;
ax3.LineWidth = 1.0;
ax3.XTickLabel = cellfun(@(x) ['\bf' x], ax3.XTickLabel, 'UniformOutput', false);