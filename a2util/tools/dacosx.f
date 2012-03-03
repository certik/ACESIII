
c This routine returns the DACOS of VALUE. If VALUE is within THRESH of 1.d0 or
c -1.d0, then VALUE is set to (-)1.d0 and the DACOS is evaluated and returned.

      double precision function dacosx(dValue,dThresh)
      implicit none
      double precision dValue, dThresh
      if (dabs(dabs(dValue)-1.d0).lt.dThresh) dValue=sign(1.d0,dValue)
      dacosx=dacos(dValue)
      return
      end

