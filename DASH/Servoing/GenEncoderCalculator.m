function encoderCalc = GenEncoderCalculator( enc_range )
%GENENCODERCALCULATOR Summary of this function goes here
%   Detailed explanation goes here

encoderCalc = @encoderCalculator;

% Pre-calculations
    center = (enc_range(2) + enc_range(1))/2;
    center_shift = (center - enc_range(1));
    new_map = @(x) (x - center_shift).*(x >= center) + (x + center_shift).*(x < center);
        

    function dist = encoderCalculator(desired, current)
        d1 = desired - current;
        d2 = new_map(desired) - new_map(current);
        
        dist = d1 .* (abs(d1) < abs(d2)) + d2 .* ~(abs(d1) < abs(d2));
    end

end

