classdef WeightedMSERegressionLayer < nnet.layer.RegressionLayer
    properties
        % You can add parameters if needed
        Epsilon
    end

    methods
       function layer = WeightedMSERegressionLayer(name, epsilon)
            % Constructor for the layer
            layer.Name = name;
            layer.Description = 'Weighted MSE Regression Layer';
            layer.Epsilon = epsilon;
        end
        function loss = forwardLoss(layer, Y, T)
            % Y: predictions
            % T: targets

            % Define weights: higher weight for larger targets
            epsilon = layer.Epsilon; % so that the weights are never zero
            weights = (T+epsilon).^2;  % penalize errors more where T is larger

            % (Optional) normalize weights to prevent exploding loss
            norm_weights = weights ./ max(weights);

            % Compute squared error
            squaredError = (Y - T).^2;

            % Apply weights
            weightedError = norm_weights .* squaredError;

            % Mean weighted loss
            loss = mean(weightedError, 'all');
        end
    end
end