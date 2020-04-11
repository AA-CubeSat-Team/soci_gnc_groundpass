function post_process(obj)

% close the log file if one exists
if (obj.logging)
    fclose(obj.log_file_id);
end

lens_s    = obj.len_passes_s;
max_len_s = max(lens_s);
min_len_s = min(lens_s);
mu_len_s  = mean(lens_s);
std_len_s = std(lens_s);
fprintf('==================================================\n');
fprintf('TOTAL passes: %d, TIME overhead: %2.2f +/- %2.2fs \n\n',...
            obj.passes,mu_len_s,std_len_s);
        
% Plot ground track afterwards if not done during sim
if (~obj.make_plot)
    figure(1)
    plot(obj.lla(2,:),obj.lla(1,:),...
        'ro','MarkerSize',3,'MarkerFaceColor','r')
end

if (obj.passes>1)
    edges = linspace(min_len_s,max_len_s,10);
    quants = [0.25,0.5,0.75];
    quants_len_s = quantile(lens_s,quants);
    lstr = cell(numel(quants),1);
    figure, hold on, grid on, box on
    histogram(lens_s,edges,'HandleVisibility','off')
    for k = 1:numel(quants)
        plot([quants_len_s(k) quants_len_s(k)],get(gca,'Ylim'),'--','LineWidth',1)
        lstr{k} = strcat(num2str(quants(k)*100),'th quantile');
    end
    legend(lstr,'Location','north','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex',...
            'FontSize',15)
    xlabel('Length of Passes [s]','FontSize',15,'Interpreter','latex')
    ylabel('Frequency','FontSize',15,'Interpreter','latex')
end

end

