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
