module thermalConstants

    implicit none

    real, parameter :: airHeatCap=1.006d3, tempAir=25.d0+273.d0, tempAir4=tempAir**4, lw=2256.d3
    real, parameter :: waterContentInit=.75, proteinContent=1.d0-waterContentInit
    real            :: skinDensityInit, QVapor

    contains

        real function airThermalCond(T, i)
        !data taken from engineeringtoolbox.com and fitted in gnuplot
        !temp input in kelvin
            implicit none

            real, intent(IN) :: T!in Kelvin
            real, parameter  :: a=-0.188521, b=0.000367259, c=0.212453
            integer :: i

            if(t < 0.)then
                print*,i,'loop #'
                error stop "Negative temp in kelvin"
            end if
            airThermalCond = a * exp(-b * (T - 273.15d0)) + c
        end function airThermalCond


        real function airDensity(T)
        !values taken from wikipedia
        !temp input in kelvin
            implicit none

            real, intent(IN) :: T !in kelvin
            real, parameter  :: pressue1ATM=101.325d3, Rspec=287.058d0

            airDensity = pressue1ATM / (Rspec * T)

        end function airDensity


        elemental real function getWaterContent(Qcurrent, watercurrent)

            implicit none

            real, intent(IN) :: Qcurrent, watercurrent

            getWaterContent = max(min(min(waterContentInit - waterContentInit * (Qcurrent / QVapor), waterContentInit), &
                                  watercurrent), 0.0)

        end function getWaterContent


        !source for three below function: Analysis of Thermal Relaxation During LaserIrradiation of Tissue, Choi & Welch 2001
        real function getSkinDensity(waterContent)
        !get dsensity of skin based upon water content. In Kg m-3

            implicit none

            real, intent(IN) :: waterContent

            getSkinDensity = 1000.d0 / (waterContent + 0.649d0*proteinContent) !new function with protein

        end function getSkinDensity


        real function getSkinHeatCap(waterContent)
        !get heat capacity for skin. In J Kg-1 K-1 

            implicit none

            real, intent(IN) :: waterContent

            getSkinHeatCap = 1000.d0*(4.2*waterContent + 1.09d0*proteinContent) !new function with protein

        end function getSkinHeatCap


        real function getSkinThermalCond(waterContent, currentDensity)
        !get thermal conductivity for skin. In W m-1 K-1 

            implicit none

            real, intent(IN) :: waterContent, currentDensity

            getSkinThermalCond = currentDensity * (6.28d-4*waterContent + 1.17d-4*proteinContent)!new function with protein

        end function getSkinThermalCond
end module thermalConstants
