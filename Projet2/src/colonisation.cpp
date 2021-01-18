#include <cstdlib>
#include <string>
#include <iostream>
#include <SDL2/SDL.h>        
#include <SDL2/SDL_image.h>
#include <fstream>
#include <omp.h>
#include <ctime>
#include <iomanip>      // std::setw
#include <chrono>
#include <mpi.h>

#include "parametres.hpp"
#include "galaxie.hpp"

int nbp=10;
 
int main(int argc, char ** argv)
{
    char commentaire[4096];
    int width, height;
    SDL_Event event;
    SDL_Window   * window;

    parametres param;


    std::ifstream fich("parametre.txt");
    fich >> width;
    fich.getline(commentaire, 4096);
    fich >> height;
    fich.getline(commentaire, 4096);
    fich >> param.apparition_civ;
    fich.getline(commentaire, 4096);
    fich >> param.disparition;
    fich.getline(commentaire, 4096);
    fich >> param.expansion;
    fich.getline(commentaire, 4096);
    fich >> param.inhabitable;
    fich.getline(commentaire, 4096);
    fich.close();

    std::cout << "Resume des parametres (proba par pas de temps): " << std::endl;
    std::cout << "\t Chance apparition civilisation techno : " << param.apparition_civ << std::endl;
    std::cout << "\t Chance disparition civilisation techno: " << param.disparition << std::endl;
    std::cout << "\t Chance expansion : " << param.expansion << std::endl;
    std::cout << "\t Chance inhabitable : " << param.inhabitable << std::endl;
    std::cout << "Proba minimale prise en compte : " << 1./RAND_MAX << std::endl;
    std::srand(std::time(nullptr));

    SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO);

    window = SDL_CreateWindow("Galaxie", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                              width, height, SDL_WINDOW_SHOWN);

    galaxie g(width, height, param.apparition_civ);
    galaxie g_next(width, height);
    galaxie_renderer gr(window);

    int deltaT = (20*52840)/width;
    std::cout << "Pas de temps : " << deltaT << " années" << std::endl;

    std::cout << std::endl;
    //parallélisation MPI
    //MPI_Init(&argc, &argv);
    
    unsigned long long temps = 0;

    std::chrono::time_point<std::chrono::system_clock> start, end1, end2;
    while (1) {
        //parallélisation MPI
        double start=MPI_Wtime();
        int rank, nbp;
        MPI_Comm globComm;
        MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm_size(MPI_COMM_WORLD,&nbp);
        const int tag = 101;
        MPI_Status status;
        if(rank==0)
        {
            gr.render(g);
            double end1=MPI_Wtime();
            double end2;
            MPI_Recv (&end2 ,1 ,MPI_DOUBLE ,1 ,tag ,globComm , &status );
            double elaps2 = end1 - start;
            double elaps1=end2-start;
            temps += deltaT;
            std::cout << "Temps passe : "
            << std::setw(10) << temps << " années"
            << std::fixed << std::setprecision(3)
            << "  " << "|  CPU(ms) : calcul " << elaps1*1000
            << "  " << "affichage " << elaps2*1000
            << "\r" << std::flush;
        }
        else
        {
            mise_a_jour(param, width, height, g.data(), g_next.data(),nbp,rank);
            g_next.swap(g);
            double end2=MPI_Wtime();
            MPI_Ssend ( &end2 , 1 , MPI_DOUBLE , 0 , tag , globComm);
        }*/
        //algo parallelisé par openMP
        /*start = std::chrono::system_clock::now();
        #pragma omp parallel num_threads(nbp)
        {
            int rank=omp_get_thread_num();
            if(rank==0)
            {
                gr.render(g);
                end1 = std::chrono::system_clock::now();
            }
            else
            {
                mise_a_jour(param, width, height, g.data(), g_next.data(),nbp,rank);
                g_next.swap(g);
                end2 = std::chrono::system_clock::now();
            }
        }
        */
        std::chrono::duration<double> elaps2 = end1 - start;
        std::chrono::duration<double> elaps1 = end2 - start;
        temps += deltaT;
        std::cout << "Temps passe : "
                  << std::setw(10) << temps << " années"
                  << std::fixed << std::setprecision(3)
                  << "  " << "|  CPU(ms) : calcul " << elaps1.count()*1000
                  << "  " << "affichage " << elaps2.count()*1000
                  << "\r" << std::flush;
        //_sleep(1000);
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
          std::cout << std::endl << "The end" << std::endl;
          break;
        }
    }
    SDL_DestroyWindow(window);
    SDL_Quit();
    //algo parallélisé par MPI
    MPI_Finalize();
    return EXIT_SUCCESS;
}
