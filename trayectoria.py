import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

import math as m

bax = 399.61879999999996
bay = -139.72506666666666
baz = -38.93926666666667
sax = 0.0005976211210047346
say = 0.0005951785069613871
saz = 0.0005879903044031379
bgx = 122.67244348991005
bgy = 13.553703191783443
bgz = -10.271588661832565
sgx = m.pi/2952
sgy = m.pi/2952
sgz = m.pi/2952
g = 9.80665

I = np.identity(3)
g_matriz = np.array([0,0,-g])

class Medicion():
	def __init__(self,ruta, archivo):
		self.__ruta = ruta		
		self.__archivo = archivo

	def time_stamp(self, fila):				# trabaja sobre toda la fila, convierte latitud a radianes
		fechahora_ts =  pd.Timestamp(fila['times'])
		return fechahora_ts

	def time(self, fila):					# trabaja sobre toda la fila, convierte el timestamp en segundos
		segundos =  ((fila['timestamp'])-self.__df["timestamp"][0])*1
		return segundos

	def prepara_dataframe(self):
		self.__df = pd.read_csv(self.__ruta+self.__archivo ,header=None)
		self.__df.columns=["timestamp","nan","nan","ax","ay","az","gx","gy","gz"]
		self.__df['timestamp']=self.__df.apply(self.time, axis=1)
		self.__df.pop('nan')
		
	def corregir_ax(self, fila):
		corregida = ((fila['ax'])- bax) * sax
		return corregida

	def corregir_ay(self, fila):
		corregida = ((fila['ay'])- bay) * say
		return corregida

	def corregir_az(self, fila):
		corregida = ((fila['az'])- baz) * saz
		return corregida

	def corregir_gx(self, fila):
		corregida = ((fila['gx'])- bgx) * sgx
		return corregida

	def corregir_gy(self, fila):
		corregida = ((fila['gy'])- bgy) * sgy
		return corregida

	def corregir_gz(self, fila):
		corregida = ((fila['gz'])- bgz) * sgz
		return corregida

	def corregir(self):
		self.__df['ax']=self.__df.apply(self.corregir_ax, axis=1)
		self.__df['ay']=self.__df.apply(self.corregir_ay, axis=1)
		self.__df['az']=self.__df.apply(self.corregir_az, axis=1)
		
		self.__df['gx']=self.__df.apply(self.corregir_gx, axis=1)
		self.__df['gy']=self.__df.apply(self.corregir_gy, axis=1)
		self.__df['gz']=self.__df.apply(self.corregir_gz, axis=1)
		return

	def alineacion(self, i_ali, f_ali):
		i_tra = f_ali+1
		self.inicio_trayecto = i_tra
		data_alin = np.array(self.__df.loc[0:f_ali,('ax','ay','az','gx','gy','gz')])
		ax_m = data_alin[:,0].mean()
		ay_m = data_alin[:,1].mean()
		az_m = data_alin[:,2].mean()
		gx_m = data_alin[:,3].mean()
		gy_m = data_alin[:,4].mean()
		gz_m = data_alin[:,5].mean()

		#tita	balanceo	roll	alabeo	gx	tita	ay	bank lateral
		#beta	elevación	pitch	cabeceo gy	beta	ax	elevation
		#gama	dirección	yaw	guiñada	gz	gama	az	heading angle
		beta = m.asin(ax_m/g)	#0.0
		#beta = m.asin(ax_m)	# Algoritmo 004

		b_s = m.sin(beta)	#ax_m
		b_c = m.cos(beta)

		tita = m.asin(-ay_m/(b_c*g))	#0.0
		#tita = m.asin(-ay_m)	#Algoritmo 004

		gama = 0.0#m.pi/4
		t_s = m.sin(tita)
		t_c = m.cos(tita)
		g_s = m.sin(gama)
		g_c = m.cos(gama)

		self.Cbn0 = np.array([[b_c*g_c, -t_c*g_s + t_s*b_s*g_c,  t_s*g_s + t_c*b_s*g_c],
				[b_c*g_s,  t_c*g_c + t_s*b_s*g_s, -t_s*g_c + t_c*b_s*g_s],
				[-b_s,		t_s*b_c,		t_c*b_c		]])

		Cnb0 = np.transpose(self.Cbn0)
		return

	def algoritmo(self):
		i_tra = self.inicio_trayecto

		dif_tiempo = self.__df['timestamp'].diff()	# calcula los delta t
		delta_t = np.array(dif_tiempo.loc[i_tra:])

		Cbn = self.Cbn0
		print(Cbn)
		Abn = np.array(self.__df.loc[i_tra:,('ax','ay','az')])
		W = np.array(self.__df.loc[i_tra:,('gx','gy','gz')])
		
		data_tray = np.zeros((len(self.__df)-i_tra,16))
		data_tray[:,(0,1,2,3,4,5,6)] = self.__df.loc[i_tra:len(self.__df),('timestamp','ax','ay','az','gx','gy','gz')]
		
		#Abn[:,0] = Abn[:,0]-ax_m
		#Abn[:,1] = Abn[:,1]-ay_m
		#Abn[:,2] = Abn[:,2]-az_m

		Dn = np.array([[0,0,0]])
		dx = 0
		dy = 0
		dz = 0

		Vn = np.array([[0,0,0]])
		Rn = np.array([[0,0,0]])

		for i in range(0, len(data_tray)):  #.index.values:
			dt=delta_t[i]

			wx = W[i,0]#-gx_m  # roll
			wy = W[i,1]#-gy_m  # pitch
			wz = W[i,2]#-gz_m  # yaw

			dx = wx*dt
			dy = wy*dt
			dz = wz*dt

			Dbn = np.array([[0,-dz,dy],
				  [dz,0,-dx],
				  [-dy,dx,0]])
			I_t = dt*I

			Cbn_k = np.dot(Cbn,I+Dbn)
			Cnb_k = np.transpose(Cbn_k)

#			Vn_k = (np.dot(Cbn_k,Abn[i]-np.dot(Cnb_k,g_matriz)))#*dt
			Vn_k = np.dot(Cbn_k,Abn[i])-np.dot(I,g_matriz)
			Vn = np.append(Vn, [dt*(Vn_k+Vn[i])], axis=0)
#			Vn = np.append(Vn, [Vn_k], axis=0) # Algoritmo 004

			Rn_k = np.dot(I,Vn[i])
			Rn = np.append(Rn, [(dt*Rn_k+Rn[i])], axis=0)

			Cbn = Cbn_k
			Cnb = np.transpose(Cbn)
			#Rn = np.delete(Rn, (0), axis=0)
			#Vn = np.delete(Vn, (0), axis=0)

		for i in range(0,len(data_tray)):
			data_tray[i,7]=Vn[i,0]
			data_tray[i,8]=Vn[i,1]
			data_tray[i,9]=Vn[i,2]
			data_tray[i,10]=Rn[i,0]
			data_tray[i,11]=Rn[i,1]
			data_tray[i,12]=Rn[i,2]
		self.df_trayectoria = pd.DataFrame(data_tray, columns=['timestamp','ax','ay','az','gx','gy','gz','vx','vy','vz','rx','ry','rz','vkx','vky','vkz'])
		return self.df_trayectoria

	def obtener_dataframe(self):
		return self.__df

	def graficar(self, var1, var2, var3, var4, var5, var6):
		x = self.df_trayectoria['timestamp'].values
#		x = self.df_trayectoria.index.values
		fig, axs = plt.subplots(3, sharex=True)
		axs[0].plot(x, self.df_trayectoria[var1].values, label=str(var1), color="blue")
		axs[0].plot(x, self.df_trayectoria[var2].values, label=str(var2), color="red")
		axs[0].legend()
		axs[0].grid()
#		axs[0].set_ylim([-20000,20000])		
		axs[0].set_title('Mediciones de acelerómetro y giróscopo')
		axs[1].plot(x, self.df_trayectoria[var3].values, label=str(var3), color="blue")
		axs[1].plot(x, self.df_trayectoria[var4].values, label=str(var4), color="red")
		axs[1].legend()
		axs[1].grid()
#		axs[1].set_ylim([-20000,20000])		
#		axs[1].set_title('Variable: '+str(var2))
		axs[2].plot(x, self.df_trayectoria[var5].values, label=str(var5), color="blue")
		axs[2].plot(x, self.df_trayectoria[var6].values, label=str(var6), color="red")
		axs[2].legend()
		axs[2].grid()
#		axs[2].set_ylim([-20000,20000])		
#		axs[2].set_title('Variable: '+str(var3))
		axs[2].set_xlabel("Tiempo [seg]")
		plt.show()
		return

	def graficar_trayecto(self, x, y):
		xtit=x
		ytit=y
		x = self.df_trayectoria[x].values
		y = self.df_trayectoria[y].values
		plt.plot(x, y, label=str(y)+" en el eje y vs. "+str(x)+" en el eje x" , color="blue")
		plt.title(ytit+" en el eje y vs. "+xtit+" en el eje x")
		plt.grid()
		plt.show()
		return

prueba = 1		
if prueba == 1:
	ruta = "./Prueba01/"
	i_ali=0
	f_ali=59000
if prueba == 2:
	ruta = "./Prueba02/"
	i_ali=0
	f_ali=6000
if prueba == 3:
	ruta = "./Prueba03/"
	i_ali=0
	f_ali=5000

archivo = "acelerometro.dat"
salida = "salida.csv"
medicion1 = Medicion(ruta, archivo)
medicion1.prepara_dataframe()
medicion1.corregir()

medicion1.alineacion(i_ali, f_ali)
df_trayectoria = medicion1.algoritmo()
print(df_trayectoria.head())

medicion1.graficar('ax','gx','ay','gy','az','gz')
medicion1.graficar('vx','rx','vy','ry','vz','rz')
medicion1.graficar_trayecto('rx', 'ry')
df = medicion1.obtener_dataframe()

pd.DataFrame(df_trayectoria).to_csv(ruta+salida, index= False)
