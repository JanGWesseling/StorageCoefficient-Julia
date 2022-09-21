module VanGenuchten

  function moistureContent(aData :: Main.Control.Types.MVGParam, ah :: Float64)
    theta = 0.0;
    try
      try
        if ah > -1.0e-4
          theta = aData.thetaS
        else
          head = ah
          # first compute |alpha * h| ** n
          help = abs(aData.alpha * head)^aData.n

          # add 1 and raise to the power m
          help = (1.0 + help) ^ aData.m

          # now compute theta
          theta = aData.thetaR + (aData.thetaS - aData.thetaR) / help
        end
      catch e
        println("?????ERROR in VanGenuchten.moistureContent: ",e)
      end
    finally
    end
    return theta
  end

  function pressureHead(aData :: Main.Control.Types.MVGParam, aTheta :: Float64)
    # calculation of h from theta
    h = 0.0
    try
      try
        help = 0.0
        if ((aTheta < aData.thetaR) || (aTheta > aData.thetaS))
          help = 1.0e3
        else
          # first calculate the inverse of the sorptivity }
          help = (aData.thetaS - aData.thetaR) / (aTheta - aData.thetaR);

          # raise to the power 1/m
          if (abs(aData.n) < 1.0e-4) || (abs(aData.m) < 1.0e-4)
            help = 1.0e28
          else
            if log10(help) / aData.m > 28.0
              help = 1.0e28
            else
              help = help ^ (1.0 / aData.m)
              if abs(help - 1.0) < 1.0e-3
                help = 1.0e28
              else
                # subtract one and raise to the power 1/n }
                if log10(help-1.0) / aData.n > 28.0
                  help = 1.0e28
                else
                  help = (help - 1.0)^(1.0 / aData.n)
                  # divide by alpha
                  help = -1.0 * Math.abs(help/aData.alpha);
                end
              end
            end
          end
        end
      catch ex
        println("ERROR in VanGenuchten.pressureHead: ",ex)
      end
    finally
    end
    return help
  end

  function conductivity(aData :: Main.Control.Types.MVGParam, ah :: Float64)
    head = ah
    term1=0.0
    term2=0.0
    alphaN=0.0
    try
      try
        if head > -0.0001
          term1 = aData.kSat
        else
          # term1 = (1 + |alfa * h| ^ n ) ^ m
          alphaH = abs(aData.alpha * head)
          logTerm = log10(alphaH)
          if aData.n * logTerm > 28.0
            term1 = 1.0e28
          else
            alphaN = alphaH^aData.n
            help = 1.0 + alphaN
            if aData.n * logTerm > 28.0
              term1 = 1.0e28
            else
              term1 = help^aData.m
            end
          end

          # term2 = |alfa * h| ^ (n-1)
          help  = aData.n - 1.0
          if help * logTerm > 28.0
            term2 = 1.0e28
          else
            term2 = alphaH^help
          end

          # the difference
          term1 = term1 - term2
          term2 = term1 * term1
          # now term2 is the nominator

          # the denominator
          term3 = 1.0 + alphaN
          if aData.m * (aData.l+2.0) * log10(term3) > 28.0
            term1 = 1.0e28;
          else
            term3 = term3 ^ (aData.m*(aData.l+2.0))
            # the conductivity in cm/d
              term1 = aData.kSat * term2 / term3;
          end
        end
        # convert to m/s
#        term1 = 0.01 * term1 / 86400.0;
      catch e
        println("?????ERROR in VanGenuchten.conductivity: ",e)
      end
    finally
    end
    return term1
  end

  function moistureCapacity(aData :: Main.Control.Types.MVGParam, ah :: Float64)
    result = 0.0
    try
      try
        if ah > -1.0e-3
          result = 0.0
        else
          head = abs(ah)
          # use analytical evaluation of capacity
          alphaH = abs(aData.alpha * head)

          # compute |alpha * h| to the power n-1
          term1 = alphaH ^  (aData.n-1.0)

          # compute |alpha*h| to the power n
          term2 = term1 * alphaH;

          # add one and raise to the power m+1
          term2 = (1.0 + term2) ^ (aData.m + 1.0)

          # divide theta-s minus theta-r by term2
          term2 = (aData.thetaS - aData.thetaR) / term2;

          # calculate the differential moisture capacity
          result = 100.0 * (aData.n * aData.m * aData.alpha * term2 * term1);
        end

      catch e
        result = 1.0e28;
        println("?????ERROR in VanGenuchten.moistureCapacity: ",e)
      end
    finally
    end
    return result;
  end

end
