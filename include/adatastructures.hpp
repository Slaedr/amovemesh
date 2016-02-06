#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#ifndef _GLIBCXX_VECTOR
#include <vector>
#endif

#ifndef _GLIBCXX_FORWARD_LIST
#include <forward_list>
#endif

#define __ADATASTRUCTURES_H

template <class T>
/// Iterator for [PList](@ref PList)
class PListIterator;

template <class T>
/// Class implementing a stack which is not compressed after arbitrary element deletion. New data can be written to unused entries rather than pushed back.
/** Uses the std::vector class as the underlying stack implementation.
 */
class PList
{
	/// Actual number of elements in the list
	int size;
	/// Number of "empty" spaces
	int ngaps;
	/// Contains data and empty spaces
	std::vector<T> data;
	/// Indices of empty spaces
	std::vector<int> freed;
	/// Indices of used elements of data
	std::forward_list<int> used;

public:
	friend class PListIterator<T>;

	PList()
	{
		ngaps = 0;
		size = 0;
	}
	
	void reserve(const int sze)
	{
		data.reserve(sze);
	}
	
	void resize(const int sze)
	{
		data.resize(sze);
		used.resize(sze);
		size = sze;
	}

	/// element access
	T& operator[](const int index)
	{
		return data[index];
	}

	/// getter
	T at(const int index) const
	{
		return data.at(index);
	}

	/// insert at end of stack
	void push_back(const T dat)
	{
		data.push_back(dat);
		used.insert_after(used.end(), data.size()-1);
		size++;
	}
	
	/// delete from end of stack
	void pop_back()
	{
		data.pop_back();
		size--;
	}
	/// delete from arbitrary location in stack
	void delete_element(const int index)
	{
		auto todelete = used.begin()+index;
		freed.push_back(*todelete);
		used.erase_after(todelete-1);
		ngaps++;
		size--;
	}
	
	/// insert at first free location (as given by freed )
	void insert_element(const T val)
	{
		if(!freed.empty())
		{
			data[freed.back()] = val;
			freed.pop_back();
			ngaps--;
		}
		else
			data.push_back(val);
	}

	/// Actually remove all deleted members
	void compress()
	{
		int j;
		for(int i = 0; i < freed.size(); i++)
		{
			data.erase(data.begin()+freed[i]);

			for(j = i+1; j < freed.size(); j++)
				if(freed[j] > freed[i])
					freed[j]--;
		}

		freed.clear();
	}

	/// output contents
	void print() const
	{
		bool toprint;
		if(freed.empty()) 
			for(int i = 0; i < data.size(); i++)
				std::cout << data[i] << " ";
		else
			for(int i = 0; i < data.size(); i++)
			{
				toprint = true;
				for(int j = 0; j < freed.size(); j++)
				{
					if(i == freed[j])
						toprint = false;
				}
				if(toprint)
					std::cout << data[i] << " ";
			}
		std::cout << endl;
	}

	/// Returns an iterator pointing to the first element
	PListIterator begin() const
	{
		// first look for the first un-freed index
		int firstpos = 0;
		for(int i = 0; i < freed.size(); i++)
			if(freed[i] == firstpos)
				firstpos++;

		PListIterator itb(*this);
	}
};

template class<T>
/// Iterator for [PList](@ref PList)
class PListIterator
{
	PList<T>& list;
	int pos;

public:
	PListIterator(PList<T>& pl, int position) : list(pl), pos(position)
	{}
};

template <class T>
struct Node
{
	T data;
	Node* next;
};

template <class T>
class CircList
{
private:
	Node<T>* start;
	int size;

public:
	CircList(T n)
	{
		start = new Node<T>;
		start->data = n;
		start->next = start;
		size = 1;
	}

	//Get the node which is "n" nodes away from start
	Node<T>* traverse(int n)
	{
		if(n < 0) { cout << "! CircList: traverse: invalid argument!\n"; return nullptr;}
		Node<T>* cur;
		cur = start;
		int i = 0;

		while(true)
		{
			if(i == n) break;
			cur = cur->next;
			i++;
		}

		return cur;
	}

	// Get last node (before start)
	Node<T>* last()
	{
		Node<T>* cur;
		cur = start;

		while(cur->next != start)
			cur = cur->next;

		return cur;
	}

	void push(T x)
	{
		size++;
		Node<T>* nnew = new Node<T>;
		nnew->data = x;
		nnew->next = start;
		last()->next = nnew;
	}

	// return pointer to first node containing data x
	Node<T>* find(T x)
	{
		if (start->data == x) return start;
		Node<T>* cur = start->next;
		while(cur != start)
		{
			if(cur->data == x) return cur;
			cur = cur->next;
		}
		cout << "! CircList: find: Data not found in list!\n";
		return nullptr;
	}

	void print_list()
	{
		Node<T>* cur;
		cur = start;
		cout << start->data << " ";

		while(cur->next != start)
		{	cur = cur->next;
			cout << cur->data << " ";
		}
		cout << endl;
	}

	~CircList()
	{
		Node<T> *cur, *nnext;
		cur = start->next;
		while(cur != start)
		{
			nnext = cur->next;
			delete cur;
			cur = nnext;
		}
		delete start;
	}
};

// This function is a cyclic permutation of consecutive integers from 'start' to 'end' (inclusive). It returns the integer (between 'start' and 'end') that is 'off' integers away from 'n' in the cyclic order.
int perm(int start, int end, int n, int off)
{
	if(n > end) { cout << "Permutation point error!\n"; return 0; }
	if(off == 0) return n;

	CircList<int> list(start);
	for(int i = start+1; i <= end; i++)
		list.push(i);

	Node<int>* nn = list.find(n);
	Node<int>* cur = nn;
	for(int i = 0; i < off; i++)
		cur = cur->next;
	return cur->data;
}

template <class T>
/// Node structure for the graph data structure
struct TNode
{
	T data;

	/// Stores pointers to neighbors
	TNode** next;
};

class AGraph
{
	TNode start;
};
