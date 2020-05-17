#include <iostream>
#include <vector>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <time.h>
#include <algorithm>
#include <string>
#include <map>

#include "polish.h"

double V_blocks=0;

std::map<std::string, Cmodule*> MapOfBlocks;


Cmodule::Cmodule(){
}

Cmodule::Cmodule(char t, int p, int left, int right){
    m_tag =t; m_parent = p;  m_left = left;  m_right = right;
    m_x = 0;    m_y = 0;    m_z = 0;
    m_w = 0;    m_d = 0;    m_h = 0;
}

Cmodule::Cmodule(int p, double w, double d, double h, std::string bname, int ptr){
    m_tag ='L'; m_parent = p;  m_left = 0;  m_right = 0;
    m_x = 0;    m_y = 0;    m_z = 0;
    m_w = w;    m_d = d;    m_h = h;
    m_name = bname; m_ptr = ptr;
}

Cmodule& Cmodule::operator = (Cmodule& b) {
  m_tag = b.m_tag;
  m_parent = b.m_parent;
  m_left = b.m_left;
  m_right = b.m_right;
  m_x = b.m_x;
  m_y = b.m_y;
  m_z = b.m_z;
  m_w = b.m_w;
  m_d = b.m_d;
  m_h = b.m_h;

  return (*this);
}

Cmodule::~Cmodule(){
}

Cmodules::Cmodules(){
}

Cmodules::~Cmodules(){

    for ( std::vector<Cmodule*>::iterator it = m_modules.begin() ; it != m_modules.end() ; it++ ) {
		delete (*it);
		}
}

void Cmodules::print () {

    std::cout<<"\n"<<"MODULES"<<"\n\n";
	for(int i= (m_modules.size()-1)/2 ; i< m_modules.size() ; i++)
	{
		m_modules[i]->print ();
	}
}

void Cmodule::print () {

        printf(" %c,  %d , %d,  %d,  %d,  %f, %f, %f, %f, %f, %f \n", m_tag, m_ptr,  m_parent, m_left, m_right, m_w, m_d, m_h, m_x, m_y, m_z);


       printf("%s \n", m_name.c_str());

        for (int i=0; i<conn_net.size(); i++)
        { printf(" net%d \n", conn_net[i]->m_num);
        }
}

void Cnets::print(){

    std::cout<<"\n"<<"NETS"<<"\n\n";
    for(int i=0;i<m_nets.size();i++)
    {
        m_nets[i]->print();
    }

}


void Cnet::print(){
    printf("net%d \n", m_num);
    for(int i=0; i<m_net.size();i++)
    {
        printf(" %s \n",m_net[i]->m_name.c_str());

    }


}


Cnet::Cnet(){
}
Cnet::~Cnet(){
}

Cnets::Cnets(){
}


Cnets::~Cnets(){
    printf ("delete nets\n");
		for ( int i = 0 ; i < m_nets.size(); i++ )
	{
		delete m_nets[i];
	}


}


int Cmodules::compute_MIV(Cnets& my_cnets){

    int total_MIV=0;
    int n1=0;
    int n2=0;
    int n3=0;
    int n4=0;


    for(int i=0;i< my_cnets.m_nets.size();i++){


        if( my_cnets.m_nets[i]->m_net.size()>1){

            int z = my_cnets.m_nets[i]->m_net[0]->m_ptr;

            int max_z = m_modules[z]->m_z;

            int min_z = m_modules[z]->m_z + m_modules[z]->m_h;


            int z1 = 0;
            int z2 = 0;



            for(int j=1; j< my_cnets.m_nets[i]->m_net.size(); j++){

                z = my_cnets.m_nets[i]->m_net[j]->m_ptr;

                z1 = m_modules[z]->m_z;
                z2 = m_modules[z]->m_z + m_modules[z]->m_h;


                if(z2<min_z){
                    min_z = z2;
                }

                if(z1>max_z){
                    max_z = z1;
                }


            }

            int miv = (max_z - min_z);

            if (miv<0){
                miv=0;
           }

            total_MIV = total_MIV + miv ;


        }

    }

    return total_MIV;






}


double Cmodules::compute_wire(Cnets& my_cnets){


    double total_length=0;

   // printf("\n\n WIRE, nets# %d \n\n", my_cnets.m_nets.size());


    for(int i=0;i< my_cnets.m_nets.size();i++){

      //      printf("\n net# %d\n\n", i);


        if( my_cnets.m_nets[i]->m_net.size()>1){

            int z = my_cnets.m_nets[i]->m_net[0]->m_ptr;

      //      printf("my_cnets.m_nets[i]->m_net[0]->m_ptr = %d \n", z-39 );

            double min_x = m_modules[z]->m_x;
            double min_y = m_modules[z]->m_y;

            double max_x = m_modules[z]->m_x + m_modules[z]->m_w;
            double max_y = m_modules[z]->m_y + m_modules[z]->m_d;

      //      printf("min_x: %g, min_y: %g, max_x: %g, max_y: %g \n ", min_x, min_y, max_x, max_y);

            /*double min_x = m_nets[i]->m_net[0]->m_x;
            double min_y = m_nets[i]->m_net[0]->m_y;

            double max_x = m_nets[i]->m_net[0]->m_x + m_nets[i]->m_net[0]->m_w;
            double max_y = m_nets[i]->m_net[0]->m_y + m_nets[i]->m_net[0]->m_d; */

            double x1 = 0;
            double x2 = 0;
            double y1 = 0;
            double y2 = 0;


            for(int j=1; j< my_cnets.m_nets[i]->m_net.size(); j++)
            {
                z = my_cnets.m_nets[i]->m_net[j]->m_ptr;

         //       printf("my_cnets.m_nets[i]->m_net[0]->m_ptr = %d \n", z-39 );

                x1 = m_modules[z]->m_x;
                x2 = m_modules[z]->m_x + m_modules[z]->m_w;
                y1 = m_modules[z]->m_y;
                y2 = m_modules[z]->m_y + m_modules[z]->m_d;

        //        printf("x1: %g, y1: %g, x2: %g, y2: %g \n ", x1, y1, x2, y2);

                if(x1<min_x){
                    min_x = x1;
                }

                if(x2>max_x){
                    max_x = x2;
                }

                if(y1<min_y){
                    min_y = y1;
                }

                if(y2>max_y){
                    max_y = y2;
                }




            }

          //  printf("min_x: %g, min_y: %g, max_x: %g, max_y: %g \n ", min_x, min_y, max_x, max_y);


            // printf("\n\n net_length: %g \n", (max_x - min_x + max_y - min_y)/2);

            total_length = total_length + (max_x - min_x + max_y - min_y)/2;


        }

    }

  //  printf("\n\n ******* total_length********: %g \n\n", total_length);


    return total_length;


}


void Cmodules::init(char* filename){

    FILE* fp;


	fp= fopen(filename,"r");

	if ( !fp )
    return ;

    char szBuf[1024];
    memset (szBuf, 0x00, sizeof(szBuf));

    if ( fgets (szBuf, sizeof(szBuf)-1, fp) == NULL ) {
        fclose (fp);
        return ;
    }

    const char* tok = " \t\r\n";

    int nr_modules = atoi (szBuf);

    printf("Number of Blocks: %d \n", nr_modules);

 //initiating supermododules

    for (int i=0; i< nr_modules-1; i++){

        char t; int p; int left; int right;

       // int num_rnd= rand() % 2;

        if (i%2 == 0){t='Y';}
        else {t='X';}



        if (i==0){p=-9;}
        else {p=i-1;}

        left=i+1;

        right=2*nr_modules-i-2;

        Cmodule* temp_module;
        temp_module = new Cmodule(t,p,left,right);

        (m_modules).push_back(temp_module);

      //  printf("%c %d %d %d\n",t, p, left, right);

    }

    //initiating modules

    for ( int i = 0 ; i < nr_modules ; i++ ) {

        memset (szBuf, 0x00, sizeof(szBuf));

        if ( fgets (szBuf, sizeof(szBuf)-1, fp) == NULL ) {

            break;
        }

        char* pch = strtok (szBuf, tok);

        if(!pch) // if nothing is in 'pch(block)' (because there is some empty line in given file)
			continue;



		std::string bname = pch;// it chould be a block name
		//printf("%s \n", pch);

		pch = strtok(NULL, tok);

		if(!pch)
			continue;

		if(strcmp(pch, "softrectangular") != 0)
			continue; // go to beginning of the 'while' loop, not go down.


		pch= strtok(NULL, tok );
		strtok(NULL, tok );
		strtok(NULL, tok );
		//pch = strtok(NULL, ")");


       // double V = atof (pch);
       // V_blocks = V_blocks + V;



   /*     int M= AR_min*10000;
        int N= AR_max*10000;
        double AR_1 = rand()% M+N ;
        double AR = AR_1/10000;
        double w = sqrt((V * AR / Layer)) ;


        double d = w /AR ; */

        double w = sqrt((atoi(pch) / init_Layer));
        double d = w;
        double h = init_Layer;
        double V = w*d*h;
        V_blocks = V_blocks + V;

        if ( (w <= 0.0) || (d <= 0.0)||(h <= 0.0) ) {

          break;

        }

        int p;

        if ((i==0)||(i==1)){p=nr_modules-2;}
        else {p=nr_modules-i-1;}

        Cmodule* temp_module;
        int ptr = nr_modules-1+i;
        temp_module = new Cmodule(p,w,d,h,bname, ptr);

        (m_modules).push_back(temp_module);

        MapOfBlocks.insert(std::make_pair(bname, temp_module));



       }   // Initialization Complete

       printf("Sum of Volumes of blocks:    %f \n", V_blocks);

}

bool open_net (Cnets& my_cnets, const char* file_name) {

    FILE* np;
    np = fopen (file_name, "r");
    if ( np == NULL ) {
    printf ("cannot open the file\n");
    return -3;
    }


    int x=1;

while(1)
	{   char buf[16];
		memset(buf, 0x00, sizeof(buf));

		if(fgets(buf, sizeof(buf)-1, np) == NULL)
		{
			break;

		}

		const char* tok = " \t\r\n";
        bool err = false;
		char* pch = strtok(buf, tok);

		if(!pch)
			continue;
        //puts(pch);


		if(strcmp(pch, "NetDegree") != 0)
			continue;

            //puts(pch);
        pch=strtok(NULL, ":" );
       // printf("%s\n",pch );

		int n_net= atoi(pch);
		//printf("%d \n",n_net);

        Cnet* temp_cnet=new Cnet;
        int n_b=0;
        Cmodule* temp_cmodule= new Cmodule;

		for(int i=0;i<n_net;i++)
        {

        //printf("%d \n", i);
		memset(buf, 0x00, sizeof(buf));

		if(fgets(buf, sizeof(buf)-1, np) == NULL)
		{
			break;

		}
        //puts(buf);
		pch=strtok(buf,tok);


        temp_cmodule=MapOfBlocks.find(pch)->second ;
        //std::cout<< temp<<"\n";


        if(temp_cmodule != NULL)
        {   n_b++;
            (temp_cnet->m_net).push_back(temp_cmodule);
            temp_cnet->m_num=x;
            (temp_cmodule->conn_net).push_back(temp_cnet);
            //std::string str = std::to_string(x);
            //temp_cnet->m_name = "net"+ std::to_string((long double)x);           //std::cout<<temp_cnet->m_net<<"\n";
        }
        else continue;

		//std:: map<std::string, CBlock*>::iterator itr;
		//itr=MapOfBlocks.find(pch);


        }


      //  std::cout<<"net# "<<x<<"    "<<temp_cnet<<"\n";

        x++;
            if(n_b > 1){
          //  printf("C\n");
                (my_cnets.m_nets).push_back(temp_cnet);

            }

        //(my_cnets.m_nets).push_back(temp_cnet);
	}

  fclose (np);

  return true;
}




void Cmodules::coordinate (int module_id){



    Cmodule* cur= m_modules[module_id];

    if ( cur->m_tag == 'L')
    return;

    m_modules[cur->m_left]->m_x = cur->m_x ;
    m_modules[cur->m_left]->m_y = cur->m_y ;
    m_modules[cur->m_left]->m_z = cur->m_z ;

    if (cur->m_tag == 'X'){

         m_modules[cur->m_right]->m_x = cur->m_x + m_modules[cur->m_left]->m_w ;
         m_modules[cur->m_right]->m_y = cur->m_y ;
         m_modules[cur->m_right]->m_z = cur->m_z ;
    }
    else if (cur->m_tag == 'Y'){

         m_modules[cur->m_right]->m_x = cur->m_x ;
         m_modules[cur->m_right]->m_y = cur->m_y + m_modules[cur->m_left]->m_d ;
         m_modules[cur->m_right]->m_z = cur->m_z ;
    }
    else if (cur->m_tag == 'Z'){

         m_modules[cur->m_right]->m_x = cur->m_x ;
         m_modules[cur->m_right]->m_y = cur->m_y ;
         m_modules[cur->m_right]->m_z = cur->m_z + m_modules[cur->m_left]->m_h ;
    }


    coordinate(cur->m_left);
    coordinate(cur->m_right);



}



void Cmodules::compute_vol (int module_id) {
  Cmodule* cur = m_modules[module_id];

  if ( cur->m_tag == 'L')
    return;

  // compute the vol of left
  compute_vol (cur->m_left);

  // compute the vol of right
  compute_vol (cur->m_right);

  // merge
  char tpl= cur->m_tag;
 // double W, D, H;

  if (tpl=='X'){
    cur->m_w = m_modules[cur->m_left]->m_w + m_modules[cur->m_right]->m_w;
    cur->m_d = std::max (m_modules[cur->m_left]->m_d, m_modules[cur->m_right]->m_d);
    cur->m_h = std::max (m_modules[cur->m_left]->m_h, m_modules[cur->m_right]->m_h);
  }
  else if (tpl=='Y'){
    cur->m_w = std::max (m_modules[cur->m_left]->m_w, m_modules[cur->m_right]->m_w);
    cur->m_d = m_modules[cur->m_left]->m_d + m_modules[cur->m_right]->m_d;
    cur->m_h = std::max (m_modules[cur->m_left]->m_h, m_modules[cur->m_right]->m_h);
  }
  else if (tpl=='Z'){
    cur->m_w = std::max (m_modules[cur->m_left]->m_w, m_modules[cur->m_right]->m_w);
    cur->m_d = std::max (m_modules[cur->m_left]->m_d, m_modules[cur->m_right]->m_d);
    cur->m_h = m_modules[cur->m_left]->m_h + m_modules[cur->m_right]->m_h;
  }

  //return W*D*H;
}





double Cmodules::volume(){
     double v=(m_modules[0]->m_w) *  (m_modules[0]->m_h) *  (m_modules[0]->m_d) ;
     return v;

}


int Cmodules::nebrmove(int b, int c){

       int p_t = m_modules[b]->m_parent;

       int p_t2 =b;

               // printf("\n err chk#2       %d\n", p_t );

        while((p_t!=0)){

//                 printf("\n err chk#3;       %d\n", p_t );

            if (p_t==c){
                    if(m_modules[p_t]->m_left == p_t2){
                        c=m_modules[p_t]->m_right;
                    }
                    else if(m_modules[p_t]->m_right == p_t2){
                        c=m_modules[p_t]->m_left;
                    }
                break;
            }



            p_t2=p_t;
            p_t= m_modules[p_t]->m_parent;
       }

    return c;

}



void Cmodules::perturb(){

    int nr_modules=m_modules.size();

  //  printf("nr_modules: %d,", nr_modules);

    int a= rand() % 2;

    // printf("\nNeighbor(1) or Rotation(0) :%d", a);

    if (a==1){  //Neighborhood Movements

       int b= rand() % (nr_modules-1) +1;
       int c= rand() % (nr_modules-1) +1;

       while (b==c){
         c= rand() % (nr_modules-1) +1;
       }



       int x = nebrmove(b,c);
       int y;

       if (x!=c){
        c=x;
       }
       else if (x==c){
            y=nebrmove(c,b);
            b=y;

       }


       int t_parent=m_modules[b]->m_parent;
       int b_p=m_modules[b]->m_parent;
       int c_p=m_modules[c]->m_parent;

       if (m_modules[b_p]->m_left == b){
            m_modules[b_p]->m_left = c;
       }
        else if (m_modules[b_p]->m_right == b){
            m_modules[b_p]->m_right = c;
       }

       if (m_modules[c_p]->m_left == c){
            m_modules[c_p]->m_left = b;
       }
        else if (m_modules[c_p]->m_right == c){
            m_modules[c_p]->m_right = b;
       }




       m_modules[b]->m_parent=m_modules[c]->m_parent;
       m_modules[c]->m_parent=t_parent;



     //  printf("\nExchange the smaller subtree with a child subtree: %d  %d", b,c);



    }


    else if (a==0){          //rotation movements


        int b= rand() % nr_modules;
        int r= rand() % 3;   // rotation type, 0 for X asis, 1 for Y axis, 2 for Z axis;
        int rnd_h= rand() % n_Layers + 1;
        double  AR_1 = rand()% 700+600 ;
        double AR = AR_1/1000;  // AR ratio from .7 to 1.3
        double V_m = m_modules[b]->m_w * m_modules[b]->m_d * m_modules[b]->m_h;



    //   printf("\n\n Rotation(0 for X axis, 1 for Y axis, @ for Z axix):%d  Module#%d", r, b);

        if (r==0){

            if(m_modules[b]->m_tag=='L'){

                m_modules[b]->m_h = rnd_h ;
                m_modules[b]->m_w = sqrt((V_m * AR /rnd_h)) ;
                m_modules[b]->m_d = m_modules[b]->m_w /AR ;



              /*  double depth= m_modules[b]->m_d;
                m_modules[b]->m_d= m_modules[b]->m_h;
                m_modules[b]->m_h= depth; */
            }
            else if(m_modules[b]->m_tag=='Y'){
                m_modules[b]->m_tag='Z';
            }
            else if(m_modules[b]->m_tag=='Z'){
                m_modules[b]->m_tag='Y';
            }
        }
        else if (r==1){

            if(m_modules[b]->m_tag=='L'){

                m_modules[b]->m_h = rnd_h ;
                m_modules[b]->m_w = sqrt((V_m * AR /rnd_h)) ;
                m_modules[b]->m_d = m_modules[b]->m_w /AR ;

               /* double width= m_modules[b]->m_w;
                m_modules[b]->m_w= m_modules[b]->m_h;
                m_modules[b]->m_h= width; */
            }
            else if(m_modules[b]->m_tag=='X'){
                m_modules[b]->m_tag='Z';
            }
            else if(m_modules[b]->m_tag=='Z'){
                m_modules[b]->m_tag='X';
            }
        }
        else if (r==2){

            if(m_modules[b]->m_tag=='L'){

                m_modules[b]->m_h = rnd_h ;
                m_modules[b]->m_w = sqrt((V_m * AR /rnd_h)) ;
                m_modules[b]->m_d = m_modules[b]->m_w /AR ;

             /* double width= m_modules[b]->m_w;
                m_modules[b]->m_w= m_modules[b]->m_d;
                m_modules[b]->m_d= width; */
            }
            else if(m_modules[b]->m_tag=='X'){
                m_modules[b]->m_tag='Y';
            }
            else if(m_modules[b]->m_tag=='Y'){
                m_modules[b]->m_tag='X';
            }
        }

        }
}



void metropolis(Cmodules& PE, Cmodules& PE_min_vol, float T, int M, Cnets& my_cnets){

    double gamma=.05*V_blocks/PE.volume();
    double nTotal = 1;
    double nAccepted = 0;
    double alpha =1;
    double beta = 20;

    for (int i=0; i<M; i++){

        Cmodules PE_temp;
        PE_temp=PE;
        PE.perturb();
        PE.compute_vol(0);
        PE.coordinate(0);

        double delta_v = alpha*(PE.volume()- PE_temp.volume()) + beta*(PE.compute_wire(my_cnets)-PE_temp.compute_wire(my_cnets));
        double delta_vmin = alpha*(PE.volume()- PE_min_vol.volume()) + beta*(PE.compute_wire(my_cnets)-PE_min_vol.compute_wire(my_cnets));


        //printf(" \n\ni=%d , %f\n",i,T);
        //PE.print();
        //double vol = PE.volume();

        float r = ((double) rand ()) /  RAND_MAX ;
        float AR_V = PE.m_modules[0]->m_w /PE.m_modules[0]->m_d;

        if((PE.m_modules[0]->m_h <= n_Layers)&&(AR_V <=AR_max)&&(AR_V>=AR_min)){

            if ( delta_v > 0 ){
	        nTotal += 1.0;}


            if ((delta_v<0)||(r< exp(-gamma* delta_v/T))){


                    if ( delta_v > 0 ){
                    nAccepted += 1.0;}

                    if ((delta_vmin<0)){

                        PE_min_vol = PE;
                    }


            }
            else {
                    PE=PE_temp;
                        }

        }
        else {
            PE=PE_temp;
        }


        //printf("\nVolume:%f \n", PE.volume());

    }

    printf ("        # acc: %g       \n", nAccepted * 100.0 / nTotal);

}


void sim_annealing(Cmodules& PE, Cmodules& PE_min_vol, Cnets& my_cnets){


    float T =100;
    float alpha=.99;
    int M = 500 *log(1+n_Layers) * PE.m_modules.size();
    float beta= .99;




    while ((T> 1)&&(M>10)){
        printf(" \n\n        M =%d        , T=%f      , minimum vol: %f , Wire_length= %f    , min_Wire_length= %f \n",
               M,T, PE_min_vol.volume(), PE.compute_wire(my_cnets), PE_min_vol.compute_wire(my_cnets));

        printf("\n MIV: %d   MIV_min: %d   \n", PE.compute_MIV(my_cnets), PE_min_vol.compute_MIV(my_cnets));
        metropolis(PE, PE_min_vol, T, M, my_cnets);
        T= alpha *T;
        M= beta*M;

    }

}

void Cmodules::plot(const char* file_name){

       if ( !file_name )
       return;

      FILE* fp = fopen (file_name, "w");

      if ( !fp )
        return;

   //   fprintf (fp, "3D BLOCKs\n");

    //  int nr=0;

     int nr_modules=m_modules.size();

     fprintf(fp, "plotcube([%g %g %g],[ %g %g %g], .1 ,[1 1 .7]);\n",   m_modules[0]->m_w,
                   m_modules[0]->m_d,  m_modules[0]->m_h,
                   m_modules[0]->m_x, m_modules[0]->m_y, m_modules[0]->m_z);


      for (int i= (nr_modules-1)/2 ; i<nr_modules; i++){

            float c[4];

            for ( int j = 0 ; j < 4 ; j++ )
            c[j] = ((double) rand ()) /  RAND_MAX ;

           // nr++;


            fprintf(fp, "plotcube([%g %g %g],[ %g %g %g], %f ,[%f %f %f]);\n",   m_modules[i]->m_w,
                   m_modules[i]->m_d,  m_modules[i]->m_h,
                   m_modules[i]->m_x, m_modules[i]->m_y, m_modules[i]->m_z, c[0], c[1], c[2], c[3]
                    );

            /*fprintf (fp, "set object %d rect from %g %g %g to %g %g %g fc rgb \"#%02x%02x%02x\" fs solid 1.0\n",
            nr, m_modules[i]->m_x, m_modules[i]->m_y, m_modules[i]->m_z,
            m_modules[i]->m_x + m_modules[i]->m_w, m_modules[i]->m_y + m_modules[i]->m_d,
            m_modules[i]->m_z + m_modules[i]->m_h,  c[0], c[1], c[2]); */
  }


 //   nr++;
//  fprintf (fp, "set object %d rect from 0,0 to %g,%g fc rgb \"white\" fs solid 1.0\n", nr, m_chip_width, m_chip_height);
//    fprintf (fp, "plot [0:%g][0:%g][0:%g] 0 title '', 0 notitle\n", m_modules[0]->m_w, m_modules[0]->m_d, m_modules[0]->m_h);
//    fprintf (fp, "pause -1 'Press any key'\n");

    fclose (fp);


}


// A = B;

Cmodules& Cmodules::operator = (Cmodules& B) {
  if ( m_modules.size() == 0 ) {
    for ( std::vector<Cmodule*>::iterator it = B.m_modules.begin() ; it != B.m_modules.end() ; it++ ) {
      Cmodule* m = new Cmodule;
      m_modules.push_back (m);
    }
  }

  for ( int i = 0 ; i < B.m_modules.size() ; i++ ) {
    *(m_modules[i]) = *(B.m_modules[i]);
  }

  return (*this);
}

int main(int argc, char* argv[])
{

	if(argc!=3)	{
		printf("usage: main<input>\n");
		return -1;
	}

	//srand (time(NULL));
	srand (0);

	Cmodules PE;
	Cnets my_cnets;

    Cmodules PE_min_vol;
   // Cmodules PE_min_vol_2;

//	PE.init(argv[1]);

	PE.init("n200.blocks");

	open_net(my_cnets, "n200.nets");


	PE.print();

	my_cnets.print();



	PE.compute_vol(0);

	PE.coordinate(0);
	PE.plot("init.plt");

	printf("\n wire: %f \n", PE.compute_wire(my_cnets));


    PE_min_vol=PE;

    sim_annealing(PE, PE_min_vol, my_cnets);



   // PE_min_vol.coordinate(0);



    printf("\n  volume: %f   \nMinimum Polish Expression:", PE.volume());




    //PE_min_vol.print();

 //   PE_min_vol.plot("final.plt");


 //    printf("\nMinimum volume: %f   \nMinimum Polish Expression:", PE_min_vol_2.volume());






    return 0;



}
