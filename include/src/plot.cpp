#include "plot.hpp"

void 
plot::draw(const std::vector<float>& x, const std::vector<double>& y){

  std::vector<float> x_plot  = x;
  std::vector<double> y_plot = y;

  assert(x_plot.size() == y_plot.size());

  for(int k = 0; k < y_plot.size()-1; k++){
    y_plot[k] = y_plot[k+1];
    x_plot[k] = x_plot[k+1];
  }

  std::vector<double>::iterator it_y = y_plot.begin();
  std::vector<float>::iterator it_x = x_plot.begin();
  for(; it_y != y_plot.end(); ){
  
    if(*it_y == 0 || *it_x > 1e7/28.0){
      it_y = y_plot.erase(it_y);
      it_x = x_plot.erase(it_x);
    }else{
      it_y++;
      it_x++;
    }
  
  }


  double y_max = y_plot[0];
  double y_min = y_plot[0];
  double x_min = 28*x_plot[0];
  double x_max = 28*x_plot[x_plot.size() - 1];
  for(int k = 0; k < y_plot.size(); k++){
    if(y_max < y_plot[k]) y_max = y_plot[k];
    if(y_min > y_plot[k] && y_plot[k] != 0.0) y_min = y_plot[k];
    y_plot[k] = std::log10(y_plot[k]);
  } 
  double delta_y = (std::log10(y_max) - std::log10(y_min))/(height);
  int delta_x = width/x_plot.size();
  if(delta_x == 0) delta_x = 1;

  std::cout.precision(2);
  for(int h = height + 2; h >= 0; h--){

    if(h == height + 1){
      std::cout << std::scientific << y_max << "|";
    }else if(h == 1){
      std::cout << std::scientific << y_min << "|";
    }else if(h == height/2+1){
      std::cout << "        |";
    }else{
      std::cout << "        |";
    }

    for(int k = 0; k < x_plot.size(); k++){

      int draw = (y_plot[k] - std::log10(y_min))/delta_y + 1;
      if( draw == h ){
        for(int l = 0; l < delta_x; l++) std::cout << "*";
      }else{
        for(int l = 0; l < delta_x; l++) std::cout << " ";
      }


    }
    std::cout << std::endl;
  }

  std::cout << "        -";
  for(int k = 0; k < x_plot.size(); k++){
    for(int l = 0; l < delta_x; l++) std::cout << "-";
  }
  std::cout << std::endl;
  std::cout << "        ";
  std::cout << std::scientific << x_min;
  int max = std::max(1.0, (double) x_plot.size() * delta_x - 14);
  for(int k = 0; k < max; k++){
    std::cout << " ";
  }
  std::cout << std::scientific << x_max << std::endl;
  std::cout << "        ";
  for(int k = 0; k < max/2+3; k++){
    std::cout << " ";
  }
  std::cout << "years ago" << std::endl;

}

