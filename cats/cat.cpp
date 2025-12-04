// g++ -o cat cat.cpp

#include <cstdio>
#include <vector>

class Cat {
  //private:
  protected:
  //public:
    int _age;
  public:
    Cat(const int age);
    virtual void Voice() const {
      printf("MeV\nI am %d\n", _age);
    }
};

Cat::Cat(const int age) {
  _age = age;
  printf("Creating cat...\n");
}

class Mainecoon : public Cat {
  public:
    int _size;
  public:
    Mainecoon(const int age) : Cat(age) {
      _size = 0;
      printf("Creating mainecoon...\n");
    }
    virtual void Voice() const {
      printf("MeV\nI am mainecoon and I am %d\n", _age);
    }
};

int main() {
  /*int a = 3;
  int* a_ptr = nullptr;
  a_ptr = new int(10);
  printf("*a_ptr = %d\n", *a_ptr);
  delete a_ptr;*/
  Cat vasyl(5);
  //vasyl.Voice();
  //printf("age = %d\n", vasyl._age);
  Mainecoon ma(25);
  //ma.Voice();
  std::vector<const Cat*> cats;
  cats.push_back(&vasyl);
  cats.push_back(&ma);
  for(const auto& it : cats) {
    it->Voice();
  }

  return 0;
}
