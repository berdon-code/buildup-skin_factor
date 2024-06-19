import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# -------------------------------------------------------------------------------------------
class BuildupTest:
    def __init__(self, q, Np, pwf0, Bo, oil_viscosity, h, porosity, compressibility, rw):
        # Operation data
        self.q = float(q)
        self.Np = float(Np)
        self.pwf0 = float(pwf0)

        # Oil data
        self.Bo = float(Bo)
        self.oil_viscosity = float(oil_viscosity)
        
        # Formation data
        self.h = float(h)
        self.porosity = float(porosity)
        self.compressibility = float(compressibility)
        self.rw = float(rw)

        # Result
        self.m = None
        self.k = None
        #self.pwf0 = None
        self.pwf_1hr = None
        

    def get_permeability(self, Horner):
        k = (-162.6 * self.q * self.oil_viscosity * self.Bo) / (Horner.m * self.h)
        k = round(k, 1)
        self.k = k
        
        self.m = Horner.m
        self.pwf_1hr = Horner.pwf_1hr
        self.reservoir_pressure = Horner.reservoir_pressure

        #print(self.pwf_1hr)

        #print(Horner.m)
        return k

    def get_skin(self):
        skin = 1.151*( (self.pwf0 - self.pwf_1hr)/self.m - np.log10(self.k/(self.porosity * self.oil_viscosity * self.compressibility * self.rw**2)) + 3.23 )
        skin = round(skin, 2)
        self.skin_pressure_drop = self.pwf_1hr - self.pwf0
        return skin


# -------------------------------------------------------------------------------------------
class HornerMethod:

    def __init__(self, pwf, t, tp, pwf0):
        
        self.pwf0 = float(pwf0)
        self.t = np.array([])
        self.pwf = np.array([])

        if self.pwf0 != 0:
            self.t = np.append(self.t, 1)
            self.pwf = np.append(self.pwf, self.pwf0)

        self.t = np.append(self.t, t)
        self.pwf = np.append(self.pwf, pwf)

        # Derivative of Pwf
        self.dPwf_dt = np.gradient(self.pwf, self.t)
        
        self.tp = tp

        self.adjusted_time = (self.tp + self.t) / self.t
        self.log_of_adjusted_time = np.log10(self.adjusted_time)

        # Derivative of Pwf
        self.dPwf = np.gradient(self.pwf, self.log_of_adjusted_time)

        #self.min_horner = 0
        self.max_horner = max(self.log_of_adjusted_time)

    #def calculate_derivative(self):
    #    dPwf = np.gradient(self.pwf, self.log_of_adjusted_time)
    #    self.dPwf = dPwf

    def interpolate_data(self):
        self.get_measured_pressure = interp1d(self.log_of_adjusted_time, self.pwf, kind="linear", fill_value="extrapolate")
        self.measured_pressures = self.get_measured_pressure(self.log_of_adjusted_time)

    def linear_fit(self, xmin, xmax):
        
        xinterval = np.linspace(xmin, xmax, 100)
        yinterval = self.get_measured_pressure(xinterval)
        
        self.extrapolate = np.polyfit(xinterval, yinterval, 1)
        
        self.log_extrapolation_points = np.linspace(0, self.max_horner, 100)
        self.extrapolated_pwf = np.polyval(self.extrapolate, self.log_extrapolation_points)
        
        self.m, self.reservoir_pressure = self.extrapolate
        self.m = round(self.m, 0)
        self.reservoir_pressure = round(self.reservoir_pressure, 0)

    def get_pwf_1hr(self, t=1):
        adjusted_t = (self.tp + t) / t
        log_of_adjusted_t = np.log10(adjusted_t)
        pwf_1hr = np.polyval(self.extrapolate, log_of_adjusted_t)

        self.pwf_1hr = pwf_1hr
        self.pwf_1hr = round(self.pwf_1hr, 0)
        self.t_1hr = log_of_adjusted_t
# -------------------------------------------------------------------------------------------
class Extrapolation_Interval:
    def __init__(self, xmin, xmax):
        self.isExtrapolated = False
        self.xmin = xmin
        self.xmax = xmax
    def extrapolate(self):
        self.isExtrapolated = True

#data = np.loadtxt("buildup2.csv", delimiter=",", skiprows=1)
#t = data[:, 0]
#pwf = data[:, 1]




#st.write(h.log_of_adjusted_time)
        

# Initializing session variables
#if "extrapolated" not in st.session_state:
#    st.session_state["extrapolated"] = False

# -------------------------------------------------------------------------------------------
if "extrapolation_interval" not in st.session_state:
    extrapolation_interval = Extrapolation_Interval(0, 1)
    st.session_state["extrapolation_interval"] = extrapolation_interval

# -------------------------------------------------------------------------------------------
st.title("Daño a la formación")
st.markdown("Con este programa se puede calcular el daño a la formación a través de los datos medidos en una prueba de Restauración de Presión utilizando el Método de Horner.")

menu = ["Datos", "Horner", "Resumen"]
choice = st.selectbox("Menu", menu)

st.markdown("---")


if choice == "Datos":

    st.markdown("#### Ingreso de datos")

    col1, col2, col3 = st.columns(3)

    with col1:
        with st.expander("Información de operación"):
            q = st.text_input("Tasa de flujo antes del período de cierre (bbl/día)", 280)
            Np = st.text_input("Producción acumulada (bbl)", 2682)
            pwf0 = st.text_input("Presión al momento de cierre (psia)", 1123)

    with col2:
        with st.expander("Información del petróleo"):
            Bo = st.text_input("Factor de volumen del petróleo (BY/BN)", 1.31)
            visc = st.text_input("Viscosidad (cp)", 2.0)

    with col3:
        with st.expander("Información de la formación"):
            h = st.text_input("Espesor neto productivo (ft)", 40)
            porosity = st.text_input("Porosidad", 0.10)
            ct = st.text_input("Compresibilidad total ($\\text{psi}^{-1}$)", 15e-6)
            rw = st.text_input("Radio de drene (ft)", 0.333)

    data_file = st.file_uploader("Subir datos", type=["csv"])

    if data_file is not None:
        data = np.loadtxt(data_file, delimiter=',', skiprows=1)

    guardar = st.button("Guardar")

    if guardar:
        st.success("Registrado correctamente")
        t = data[:, 0]
        pwf = data[:, 1] 

        tp = (float(Np)/float(q))*24
        tp = round(tp, 0)

        
        st.session_state["buildup_test"] = BuildupTest(q, Np, pwf0, Bo, visc, h, porosity, ct, rw)

        st.session_state["horner"] = HornerMethod(pwf, t, tp, pwf0)

        print(tp)
        

elif choice == "Horner":
    st.markdown("#### Método de Horner")
    

    if "horner" in st.session_state:
    
        buildup_test = st.session_state["buildup_test"]
        horner = st.session_state["horner"]

        fig, ax = plt.subplots()

        #ax.plot(horner.t, horner.pwf, 'o--') 
        #ax.plot(horner.t, horner.dPwf_dt, 'o--')

        ax.plot(horner.log_of_adjusted_time, horner.pwf, 'o--', label="Data")
        #ax.plot(horner.log_of_adjusted_time[::-1], horner.pwf, 'o--') #inverted
        #ax.plot(horner.log_of_adjusted_time, horner.pwf, 'o--') #inverted

        #ax.plot(horner.log_of_adjusted_time, horner.dPwf, 'o--', label="Derivada")

        ax.set_xlim([0, horner.max_horner + 0.1])

        ax.set_xlabel("$\log(\\frac{(%g+t)}{t})$"%(horner.tp))
        #ax.set_xlabel("(tp+t)/t")
        ax.set_ylabel("Presión de fondo (psi)")
        
        ax.invert_xaxis()
        ax.grid()
        ax.legend()
        st.pyplot(fig)

        st.markdown("---")
        # Call extrapolation interval
        extrapolation_interval = st.session_state["extrapolation_interval"]
        if not extrapolation_interval.isExtrapolated:
            extrapolation_interval.xmax = round(horner.max_horner, 2)

        st.markdown("#### Linear fit")
        col1, col2 = st.columns(2)
        with col1:
            xmin = st.text_input("Inicio del intervalo", extrapolation_interval.xmin)
        with col2:    
            xmax = st.text_input("Final del intervalo", extrapolation_interval.xmax)
        fit = st.button("Extrapolar") 
        
        if fit:
            extrapolation_interval = Extrapolation_Interval(xmin, xmax)
            extrapolation_interval.extrapolate()
            st.session_state["extrapolation_interval"] = extrapolation_interval

        if st.session_state["extrapolation_interval"].isExtrapolated:
            horner.interpolate_data()  
            horner.linear_fit(float(xmin), float(xmax)) 

            horner.get_pwf_1hr()
            
            fig, ax = plt.subplots()
            horner = st.session_state["horner"]
            ax.plot(horner.log_of_adjusted_time, horner.pwf, 'o--')
            ax.plot(horner.log_extrapolation_points, horner.extrapolated_pwf, '-', label="$p_{wf}=%g[\\log(\\frac{(%g+t)}{t})] + %g$"%(horner.m, horner.tp, horner.reservoir_pressure))
            
            # Plot reservoir pressure
            ax.plot(0, horner.reservoir_pressure, 'o', label=f"Presión del yacimiento: {horner.reservoir_pressure}")

            # Plot extrapolated pwf_1hr
            ax.plot(horner.t_1hr, horner.pwf_1hr, 'o', label=f"Presión extrapolada a 1 hora: {horner.pwf_1hr}")

            # Plot extrapolated pwf_1hr
            ax.plot(horner.t_1hr, horner.pwf0, 'o', label=f"Pwf($\Delta t = 0$): {horner.pwf0}")

            a = horner.t_1hr, horner.t_1hr
            b = horner.pwf0, horner.pwf_1hr

            print(a)
            print(b)

            ax.plot(a, b, '--', label="Caída de presión por daño a la formación")
            #ax.plot(a, b, '-', label="Skin")

            ax.set_xlabel("$\log(\\frac{(%g+t)}{t})$"%(horner.tp))
            ax.set_ylabel("Presión de fondo (psia)")

            ax.invert_xaxis()
            
            ax.grid()
            ax.legend()
            st.pyplot(fig)

            st.markdown("---")
            st.markdown("#### Permeabilidad")
            k = buildup_test.get_permeability(horner)
            st.markdown("##### $ k = \\frac{B_oq\\mu_o}{mh} = \\frac{(%s)(%s)(%s)}{(%s)(%s)} = %s \\text{ mD} $"%(buildup_test.Bo, buildup_test.q, buildup_test.oil_viscosity, buildup_test.m, buildup_test.h, k))

            st.markdown("---")
            st.markdown("#### Factor de daño")
            s = buildup_test.get_skin()
            st.markdown("##### $S = 1.151[\\frac{p_{wf}(t_p) - p_{1h}}{m} - \\log{\\frac{k}{\\Phi c_tr^2_w}} + 3.23 ] $")
            st.markdown("##### $S = 1.151[\\frac{(%s) - %s}{%s} - \\log{\\frac{%s}{(%s)(%s)(%s)^2}} + 3.23] = %s $" %(buildup_test.pwf0, buildup_test.pwf_1hr, buildup_test.m, buildup_test.k, buildup_test.porosity, buildup_test.compressibility, buildup_test.rw, s))


    else:
        st.warning("Ingresa los datos primero")

# elif choice == "Permeabilidad":
#     st.markdown("#### Permeabilidad")

#     buildup_test = st.session_state["buildup_test"]
#     horner = st.session_state["horner"]

#     col1, col2 = st.columns(2)
#     with col1:
#         k = st.text_input("Permeabilidad (mD)", buildup_test.get_permeability(horner))


# elif choice == "Skin":
#     st.markdown("#### Factor de daño")
    
#     buildup_test = st.session_state["buildup_test"]

#     horner = st.session_state["horner"]

#     col1, col2 = st.columns(2)
#     with col1:
#         s = st.text_input("Factor de daño", buildup_test.get_skin())
#         #s = st.text_input("Factor de daño", "Daño")

elif choice == "Resumen":
    st.markdown("#### Resumen")

    if st.session_state["extrapolation_interval"].isExtrapolated:

        buildup_test = st.session_state["buildup_test"]
        horner = st.session_state["horner"]

        col1, col2, col3 = st.columns(3)
        with col1:
            reservoir_pressure = st.text_input("Presión del yacimiento (psi)", buildup_test.reservoir_pressure)
            k = st.text_input("Permeabilidad (mD)", buildup_test.get_permeability(horner))

        with col2:
            pwf_1hr = st.text_input("Presión extrapolada a una hora (psi)", buildup_test.pwf_1hr)
            s = st.text_input("Daño", buildup_test.get_skin())

        with col3:
            pwf0 = st.text_input("Presión medida al tiempo de cierre (psi)", buildup_test.pwf0)
            skin_pressure_drop = st.text_input("Caída de presión por daño a la formación", buildup_test.skin_pressure_drop)
    
    else:
        st.warning("Ingresa los datos primero")







