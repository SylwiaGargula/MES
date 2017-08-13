#include <iostream>
#include <cstdlib>

using namespace std;


struct wezel
{
	double r0; //wspolrzedna
	int stan; //0 brak,  lub tylko 2 koniec 
	int ID; //indeks elementu w siatce
	int waga; //waga wezla =1
	
};

struct element
{
	 
	double k; //wspolczynnik
	double alfa; //warunek 2
	double tot; //temperatura otoczenia

	double c; //cieplo wlasciwe
	double deltar; //skok r
	double rmax; //maksymalny promien
	double ro; //gestosc
	double deltatau;
	double temperaturapoczatkowa1;
	double temperaturapoczatkowa2;

	wezel wezly[2]; //wezly dwa
	double H[2][2]; //macierz pojemnosci cieplnej
	double P[2]; // macierz obciazen
	
	void oblicz_macierze_lokalne()
	{
		//macierz H
		for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			H[i][j] = 0;

		double ksi = 0.5573;
		//funkcje ksztaltu
		double Ni[2];
		double Nj[2];

		Ni[0] = 0.5*(1 - ksi);
		Ni[1] = 0.5*(1 - (-ksi));

		Nj[0]= 0.5*(1 + ksi);
		Nj[1] = 0.5*(1 + (-ksi));

		//r
		double r[2];
		r[0] = (Ni[0] * wezly[0].r0) + (Nj[0] * wezly[1].r0);
		r[1] = (Ni[1] * wezly[0].r0) + (Nj[1] * wezly[1].r0);
		

		H[0][0] = ((k / deltar)*((r[0] * wezly[0].waga) + (r[1] * wezly[1].waga))) + ((c*ro*deltar) / deltatau)*((Ni[0] * Ni[0] * r[0] * wezly[0].waga) + (Ni[1] * Ni[1] * r[1] * wezly[1].waga));
		H[0][1] = ((-(k / deltar))*((r[0] * wezly[0].waga) + (r[1] * wezly[1].waga))) + ((c*ro*deltar) / deltatau)*((Ni[0] * Nj[0] * r[0] * wezly[0].waga) + (Ni[1] * Nj[1] * r[1] * wezly[1].waga));
		H[1][0] = ((-(k / deltar))*((r[0] * wezly[0].waga) + (r[1] * wezly[1].waga))) + ((c*ro*deltar) / deltatau)*((Ni[0] * Nj[0] * r[0] * wezly[0].waga) + (Ni[1] * Nj[1] * r[1] * wezly[1].waga));
		
		if( wezly[1].stan==2)
		H[1][1] = ((k / deltar)*((r[0] * wezly[0].waga) + (r[1] * wezly[1].waga))) + ((c*ro*deltar) / deltatau)*((Nj[0] * Nj[0] * r[0] * wezly[0].waga) + (Nj[1] * Nj[1] * r[1] * wezly[1].waga))+2*alfa*rmax;
		else
		H[1][1] = ((k / deltar)*((r[0] * wezly[0].waga) + (r[1] * wezly[1].waga))) + ((c*ro*deltar) / deltatau)*((Nj[0] * Nj[0] * r[0] * wezly[0].waga) + (Nj[1] * Nj[1] * r[1] * wezly[1].waga));
		
		//macierz P
		for (int i = 0; i < 2; i++)
			P[i] = 0;
		
		P[0] = ((-(c*ro*deltar)) / deltatau)*(((Ni[0] * temperaturapoczatkowa1 + Nj[0] * temperaturapoczatkowa2)*Ni[0] * r[0] * wezly[0].waga) + ((Ni[1] * temperaturapoczatkowa1 + Nj[1] * temperaturapoczatkowa2)*Ni[1] * r[1] * wezly[1].waga));
		if (wezly[1].stan == 2)
			P[1] = ((-(c*ro*deltar)) / deltatau)*(((Ni[0] *temperaturapoczatkowa1 + Nj[0] *temperaturapoczatkowa2)*Nj[0] * r[0] * wezly[0].waga) + ((Ni[1] *temperaturapoczatkowa1 + Nj[1] * temperaturapoczatkowa2)*Nj[1] * r[1] * wezly[1].waga)) - 2 * alfa*rmax*tot;
		else
			P[1] = ((-(c*ro*deltar)) / deltatau)*(((Ni[0] *temperaturapoczatkowa1 + Nj[0] *temperaturapoczatkowa2)*Nj[0] * r[0] * wezly[0].waga) + ((Ni[1] * temperaturapoczatkowa1 + Nj[1] *temperaturapoczatkowa2)*Nj[1] * r[1] * wezly[1].waga));
		
	

	}
};

struct siatka
{
	int ne; //liczba elementow
	int nh; //liczba wierzcholkow
	double **GH; //globalna 
	double *GP;
	double *Gt; //wynik
	element *elementy;
	

	void stworz_macierze()
	{

		//macierz GH
		GH = new double*[nh];
		for (int i = 0; i < nh; i++)
		{
			GH[i] = new double[nh];
		}

		for (int i = 0; i < nh; i++)
		{
			for (int j = 0; j < nh; j++)
			{
				GH[i][j] = 0;
			}
		}

		//macierz GP Gt
		GP = new double[nh];
		Gt = new double[nh];
		for (int i = 0; i < nh; i++)
		{
			GP[i] = 0;
			Gt[i] = 0;
		}

		//sumowanie macierz
		//GH
		for (int i = 0; i < ne; i++)
		{//ustalanie miejsca w duzej macierzy
			GH[elementy[i].wezly[0].ID - 1][elementy[i].wezly[0].ID - 1] += elementy[i].H[0][0];
			GH[elementy[i].wezly[1].ID - 1][elementy[i].wezly[1].ID - 1] += elementy[i].H[1][1];
			GH[elementy[i].wezly[1].ID - 1][elementy[i].wezly[0].ID - 1] += elementy[i].H[1][0];
			GH[elementy[i].wezly[0].ID - 1][elementy[i].wezly[1].ID - 1] += elementy[i].H[0][1];
	
		}
	
		//GP
		for (int i = 0; i < ne; i++)
		{
			GP[elementy[i].wezly[0].ID - 1] += elementy[i].P[0];
			GP[elementy[i].wezly[1].ID - 1] += elementy[i].P[1];
		}
		cout.precision(15);
	//	cout << "Macierz GH" << endl;
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0; j < nh; j++)
			{
		//		cout << GH[i][j] << " ";
			}
		//	cout << endl;
		}
		cout << endl;
	//	cout << "Macierz GP" << endl;
	//	for (int i = 0; i < nh; i++)
	//	{
	//	cout << GP[i] << " ";
		
		//	cout << endl;
	//	}
		cout << endl;
	}

	void gauss()
	{
		for (int i = 0; i < nh - 1; i++)
		{
			for (int j = i + 1; j < nh; j++)
			{

				double q = -GH[j][i] / GH[i][i];
				for (int k = i; k <= nh; k++)
				{
					GH[j][k] += q* GH[i][k];
				}
				GP[j] += q*GP[i];
			}
		}

		Gt[nh] = GP[nh - 1] / GH[nh - 1][nh - 1];
		for (int i = nh - 1; i >= 0; i--)
		{
			double pom = 0;
			for (int j = i + 1; j<nh; j++)
			{
				pom += GH[i][j] * Gt[j];
			}
			Gt[i] = (GP[i] - pom) / GH[i][i];

		}
		for (int i = 0; i < nh; i++)
		{
			Gt[i] *= -1;
		}
		cout << "Macierz Gt" << endl;
		for (int i = 0; i < nh; i++)
		{
			cout << "Gt " << i + 1 << " = " << Gt[i] << endl;
		}
		
	}
	
	};
	


	int main()
	{
		
		wezel * listawezlow = new wezel[5];
	for (int i = 0; i < 5; i++)
	{
		listawezlow[i].r0 = 0.02*i;
		listawezlow[i].stan = 0;
		listawezlow[i].ID = i+1;
		listawezlow[i].waga = 1;
		

	}
	listawezlow[4].stan = 2;
	
	

	element *listaelementow = new element[4];
	
	for (int i = 0; i < 4; i++)
	{
		listaelementow[i].wezly[0] = listawezlow[i];
		listaelementow[i].wezly[1] = listawezlow[i+1];
		listaelementow[i].rmax = 0.08;
		listaelementow[i].deltar = 0.02;
		listaelementow[i].alfa = 7;
		listaelementow[i].deltatau = 50;
		listaelementow[i].tot = 253;
		listaelementow[i].temperaturapoczatkowa1 = 293;
		listaelementow[i].temperaturapoczatkowa2 = 293;
		
	}
	//ziemia
	for (int i = 0; i < 2; i++)
	{
		listaelementow[i].ro = 1800;
		listaelementow[i].c = 1200;
		listaelementow[i].k = 0.9;
	}
	//izolacja
	listaelementow[2].ro = 250;
	listaelementow[2].c = 1460;
	listaelementow[2].k = 0.07;

	//doniczka
	listaelementow[3].ro = 1800;
	listaelementow[3].c = 1200;
	listaelementow[3].k = 0.85;

	for (int i = 0; i < 4;i++)
	listaelementow[i].oblicz_macierze_lokalne();


	siatka Siatka;
	Siatka.ne = 4;
	Siatka.nh = 5;
	Siatka.elementy = listaelementow;
	Siatka.stworz_macierze();
	Siatka.gauss();

	for (int j = 0; j < 400;j++)
	{
		for (int i = 0; i < 4; i++)
		{
			Siatka.elementy[i].temperaturapoczatkowa1 = Siatka.Gt[Siatka.elementy[i].wezly[0].ID - 1];
			Siatka.elementy[i].temperaturapoczatkowa2 = Siatka.Gt[Siatka.elementy[i].wezly[1].ID - 1];
			Siatka.elementy[i].oblicz_macierze_lokalne();
		}

		Siatka.stworz_macierze();
		Siatka.gauss();
	}

		system("PAUSE");
		return 0;
	}
	