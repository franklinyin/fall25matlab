function plot_scope(out, scopeNames, xlabel_str, ylabel_str, title_str, filename)
    % Unified function to plot one or multiple scopes
    % scopeNames can be a cell array (even with one element) or a string
    
    % Convert string to cell array if needed
    if ischar(scopeNames) || isstring(scopeNames)
        scopeNames = {scopeNames};
    end
    
    figure;
    hold on; grid on;
    for i = 1:numel(scopeNames)
        name = scopeNames{i};
        fromSim = out.get(name);
        for j = 1:fromSim.numElements
            toPlot = fromSim.get(j);
            % If only one scope, use simpler display name
            if numel(scopeNames) == 1
                plot(toPlot.Values.Time, toPlot.Values.Data, 'DisplayName', name);
            else
                plot(toPlot.Values.Time, toPlot.Values.Data, 'DisplayName', [name]);
            end
        end
    end
    xlabel(xlabel_str); ylabel(ylabel_str);
    title(title_str);
    legend('show');
    saveas(gcf, filename);
end

