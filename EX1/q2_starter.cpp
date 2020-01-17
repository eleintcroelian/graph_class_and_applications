#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>


//alias pointer for convenience
using book_ptr = std::vector<std::vector<string>>*

//Helper function to split string into a vector based on a splitting character.
//This function is used to load the book.
const std::vector<std::string> split(const std::string& s, const char& c)
{
	std::string buff{""};
	std::vector<std::string> v;
	
	//iterate character by character and push the buffer to the vector if the splitting char is reached.
	for(auto n:s)
	{
		if(n != c) buff+=n; else
		if(n == c && buff != "") { v.push_back(buff); buff = ""; }
	}
	if(buff != "") v.push_back(buff);
	
	return v;
}

//This class represents a book as a vector of vector of strings,
//where each vector of strings is sentence in the book.
class Book:
{
	std::vector<std::string> temp_book_;
    std::vector<std::vector<std::string>> book_;
	
	public:
	//Constructor for the book class.
	//Takes in a .txt file and splits it into a
	//std::vector<std::vector<string>>
	Book(std::string filename) 
	{
		//Read in the book as a vector of strings
		//where each string is a sentence.
		std::string sent;
		std::ifstream book_file(filename);
		std::vector<std::string> temp;

		//Read the file line by line.
		while(std::getline(book_file, sent))
		{
			temp_book_.push_back(sent);
		}
		book_file.close();
		
		//Split each sentence by ' '
		for (unsigned int i = 0; i< temp_book_.size(); i++)
		{
			//Skip empty lines
			if (temp_book_[i].length()==0){
				continue;
			}
			else{
				temp = split(temp_book_[i], ' ');
				book_.push_back(temp);
			}
			
		}
		
	}
	

	//return an iterator pointing at the start of the book
	BookIter begin()
	{
		//Your code here
	} 
	
	//return an iterator pointing at the end of the book
	BookIter end()
	{
		//Your code here	
	}

	private:
		std::vector<std::vector<string>> book_;
	
}

class BookIter:
{
	public:
	

	//Increments to the next word in the book class.
	BookIter operator++()
	{
		//Your code here
	}
	
	//Defines equality between two iterators
	bool operator==(BookIter book_iter)
	{
		//Your code here		
	}
	
	//Defines inequality between two iterators
	bool operator!=(BookIter book_iter)
	{
		//Your code here		
	}
	
	//Dereference operator
	std::string operator*()
	{
		//Your code here			
	}
	
	
	private:
		friend class Book;
		//Your code here
	
		//Private constructor that can be accessed by the Book class.
		BookIter(){
			//Your code here
		}
}

int main()
{
	//Read in the book
	Book moby_dick("moby_dick.txt");

	//find the longest word in the book.
	//Your code here
	std::cout<< longest_word << std::endl;
    return 0;
}