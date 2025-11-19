function plot_scopes_separate(out, scopeNames, xlabel_str, ylabel_str, title_prefix, filename_prefix)
    for i = 1:numel(scopeNames)
        name = scopeNames{i};
        fromSim = out.get(name);
        
        figure;
        hold on; grid on;
        for j = 1:fromSim.numElements
            toPlot = fromSim.get(j);
            plot(toPlot.Values.Time, toPlot.Values.Data, 'DisplayName', toPlot.Name);
        end
        xlabel(xlabel_str); ylabel(ylabel_str);
        title([title_prefix ' - ' name]);
        legend('show');
        saveas(gcf, [filename_prefix '_' name '.png']);
    end
end

