#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string>
#include "../include/carpetas.h"

int ExisteFichero(char *nombreFichero)
{
	/*Función que verifica si existe o no el fichero llamado <nombreFichero> */

	struct stat stFileInfo;
	int intStat;

  	// Llamada para recoger los atributos del fichero
	intStat = stat(nombreFichero,&stFileInfo);
	if(intStat == 0)
		return 1;
	else
		return 0;
}

int CrearCarpeta(int LIMITE_FICHEROS)//, char _inputname)
{
	/*Esta función creará una carpeta con todos los permisos abiertos (0777), respetando la siguiente nomenclatura: <nombregenerico>#, es decir, si nombregenerico="PIC", creará: PIC1, PIC2, PIC3.... en orden, comprobando cual es la última versión.
	Posteriormente, el proceso se moverá a esa carpeta.
	Valores que se le pasan:
		LIMITE_FICHEROS: numero de ficheros limite que admitiremos
	Valores de retorno:
		1: se ha creado correctamente;
		0: no se ha podido crear porque se excedio el numero de ficheros_limite
	*/

	char foldername[40]="PIC";//_inputname;	//Nombre genérico que asociamos a la carpeta #: <nombregenerico>#
	char nombrecarpeta[40];		//Nombre completo de la carpeta
	int k;

		
	k=1;
	sprintf(nombrecarpeta, "%s%d", foldername, k);
	
	while(ExisteFichero(nombrecarpeta) && k<=LIMITE_FICHEROS)
	{
		sprintf(nombrecarpeta, "%s%d", foldername, k);
		k++;
	}
	
	if(k<=LIMITE_FICHEROS)
	{
		if(mkdir(nombrecarpeta, 0777)!=0)
		{
			printf ("mkdir() failure; terminating");
 			exit(1);
		}
		//printf("Hemos creado la carpeta con el nombre: %s\n", nombrecarpeta);
		//Nos movemos al directorio que hemos creado
		if(chdir(nombrecarpeta))
		{	
			printf ("chdir() failure; terminating");
 			exit(1);
		}
		//printf ("%s\n",get_current_dir_name());
		return 1;
	}
	else
	{	printf("No hemos podido crear la carpeta al superar el limite de carpetas\n");
		return 0;
	}
}

