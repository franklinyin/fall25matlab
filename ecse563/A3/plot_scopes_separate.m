function plot_scopes_separate(out, scopeNames, xlabel_str, ylabel_str, title_prefix, filename_prefix)
    % Plot multiple scopes side-by-side in the same figure
    % scopeNames should be a cell array of scope names
    
    figure;
    numScopes = numel(scopeNames);
    
    for i = 1:numScopes
        name = scopeNames{i};
        fromSim = out.get(name);
        
        % Create subplot (left and right)
        subplot(1, numScopes, i);
        hold on; grid on;
        
        for j = 1:fromSim.numElements
            toPlot = fromSim.get(j);
            plot(toPlot.Values.Time, toPlot.Values.Data, 'DisplayName', name);
        end
        
        xlabel(xlabel_str); 
        ylabel(ylabel_str);
        title([title_prefix ' - ' name]);
        legend('show');
    end
    
    saveas(gcf, [filename_prefix '.png']);
end

