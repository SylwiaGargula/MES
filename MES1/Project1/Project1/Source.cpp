#include <iostream>
#include <cstdlib>

using namespace std;


struct wezel
{
	double x; //wspolrzedna
	int stan; //0 brak, 1 poczatek grzany, 2 koniec 
	int ID; //indeks elementu w siatce
};

struct element
{
	 
	double L; //dlugosc odcinka
	double S; //pole
	double k; //wspolczynnik
	double alfa; //warunek 2
	double q;// strumien warunek 1 
	double tot; //temperatura otoczenia
	wezel wezly[2]; //wezly dwa
	double H[2][2]; //macierz pojemnosci cieplnej
	double P[2]; // macierz obciazen
	
	void oblicz_macierze_lokalne()
	{
		//macierz H
		double C = (S*k) / L;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				H[i][j] = C;
				if (i != j)
				{
					H[i][j] *= -1;
				}
			}
			
		}

		//macierz P
		for (int i = 0; i < 2; i++)
			P[i] = 0;

		for (int i = 0; i < 2; i++)
		{
			if (wezly[i].stan == 1)
			{
				P[0] = q*S;
			}
			else if (wezly[i].stan == 2)
			{
				P[1] = (-alfa)*tot*S;
				H[1][1] += alfa*S;
			}
		}
		

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
		cout << "Macierz GH" << endl;
		for (int i = 0; i < nh; i++)
		{
			for (int j = 0; j < nh; j++)
			{
				cout << GH[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "Macierz GP" << endl;
		for (int i = 0; i < nh; i++)
		{
		cout << GP[i] << " ";
		
			cout << endl;
		}
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
		
		int ilosc_elementow = 0;
		int ilosc_wierzcholkow = 0;

		cout << "Wyznaczanie ustalonego pola temperatury w precie" << endl;
		cout << "Podaj ilosc elementow" << endl;
		cin >> ilosc_elementow;
		ilosc_wierzcholkow = ilosc_elementow + 1;
		
		
		element * elementyG = new element[ilosc_elementow];
		wezel * wezlyG = new wezel[ilosc_wierzcholkow];

		cout << "Podaj dla kazdego wierzcholka wspolrzedna x, stan i ID" << endl;
		for (int i = 0; i < ilosc_wierzcholkow; i++)
		{
			cin >> wezlyG[i].x;
			cin >> wezlyG[i].stan;
			cin >> wezlyG[i].ID;
		}
		cout << endl;

		cout << "Podaj dla kazdego elementu dlugosc odcinka L, pole S, wspolczynnik k, alfa(jesli nie dotyczy to 1, warunek 2(koniec)), strumien q(jesli nie dotyczy 1, warunek 1(grzanie)), temperatura otoczenia too i numery jego wezlow" << endl;
		for (int i = 0; i < ilosc_elementow; i++)
		{
			int id_wezel1 = 0;
			int id_wezel2 = 0;
			cin >> elementyG[i].L;
			cin >> elementyG[i].S;
			cin >> elementyG[i].k;
			cin >> elementyG[i].alfa;
			cin >> elementyG[i].q;
			cin >> elementyG[i].tot;
			cin >> id_wezel1;
			cin >> id_wezel2;
			for (int j = 0; j < ilosc_wierzcholkow; j++)
			{
				if (wezlyG[j].ID == id_wezel1)
				{
					elementyG[i].wezly[0] = wezlyG[j];
				}
				else if(wezlyG[j].ID == id_wezel2)
				{
					elementyG[i].wezly[1] = wezlyG[j];
				}
			}
			elementyG[i].oblicz_macierze_lokalne();
		}
		cout << endl;
		siatka siatka = { ilosc_elementow,ilosc_wierzcholkow };
		
		siatka.elementy = elementyG;
		siatka.stworz_macierze();
		siatka.gauss();
		
		
		system("PAUSE");
		return 0;
	}
	