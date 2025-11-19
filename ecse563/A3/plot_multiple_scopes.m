function plot_multiple_scopes(out, scopeNames, xlabel_str, ylabel_str, title_str, filename, useInterpreter)
    figure;
    hold on; grid on;
    for i = 1:numel(scopeNames)
        name = scopeNames{i};
        fromSim = out.get(name);
        for j = 1:fromSim.numElements
            toPlot = fromSim.get(j);
            plot(toPlot.Values.Time, toPlot.Values.Data, 'DisplayName', [name ' - ' toPlot.Name]);
        end
    end
    xlabel(xlabel_str); ylabel(ylabel_str);
    title(title_str); saveas(gcf, filename);
    if useInterpreter
        legend('show', 'Interpreter', 'none');
    else
        legend('show');
    end
end

