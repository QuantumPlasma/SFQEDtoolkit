#pragma once

template <class myClass> class Vector3D {
private:
	myClass x;
	myClass y;
	myClass z;
	std::string vectorName;

public:
	//constructors
	//default constructor
	Vector3D() {

	}

	Vector3D(myClass defaultVal) {
		this->x = defaultVal;
		this->y = defaultVal;
		this->z = defaultVal;
	}
	
	Vector3D(myClass defaultVal, std::string name) : Vector3D(defaultVal){
		this->vectorName = name;
	}
	
	Vector3D(myClass x, myClass y, myClass z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	
	Vector3D(myClass x, myClass y, myClass z, std::string name) : Vector3D(x, y, z) {
		this->vectorName = name;
	}

	operator double*(){
		return new double[3]{this->x, this->y, this->z};
	}

	//setters
	void setX(myClass x) {
		this->x = x;
	}
	void setY(myClass y) {
		this->y = y;
	}
	void setZ(myClass z) {
		this->z = z;
	}
	void setName(std::string name) {
		this->vectorName = name;
	}

	//getters
	myClass getX() const {
		return x;
	}
	myClass getY() const {
		return y;
	}
	myClass getZ() const {
		return z;
	}
	myClass getSqModule() const {
		return *this * *this;
	}
	std::string getName() {
		return vectorName;
	}

	//overload of the * operator to make it look like the dot product
	myClass operator*(const Vector3D& vec) const {
		return (this->x * vec.x) + (this->y * vec.y) + (this->z * vec.z);
	}

	//overload of the * operator to make it look like the multiplication by scalar
	Vector3D operator*(myClass scalar) const {
		return Vector3D(this->x * scalar, this->y * scalar, this->z * scalar);
	}

	//overload of the / operator to make it look like the cross product
	Vector3D operator/(const Vector3D& vec) const {
		return Vector3D(this->y * vec.z - this->z * vec.y, this->z * vec.x - this->x * vec.z, this->x * vec.y - this->y * vec.x);
	}

	//overload of the + operator to make it look like a proper vector summation operator
	Vector3D operator+(const Vector3D& vec) const {
		return Vector3D(this->x + vec.x, this->y + vec.y, this->z + vec.z);
	}

	Vector3D operator-(const Vector3D& vec) const {
		return Vector3D(this->x - vec.x, this->y - vec.y, this->z - vec.z);
	}

	//overload of the += operator to make it look like a proper vector summation operator
	Vector3D operator+=(const Vector3D& vec) {
		this->x += vec.x;
		this->y += vec.y;
		this->z += vec.z;

		return *this;
	}

	Vector3D operator-=(const Vector3D& vec) {
		this->x -= vec.x;
		this->y -= vec.y;
		this->z -= vec.z;

		return *this;
	}

	//overloading ostream operator
	friend std::ostream& operator<<(std::ostream& os, const Vector3D& vec) {
		std::string tmp = vec.vectorName.empty() ? "Unnamed" : vec.vectorName;
		return os << tmp << " ==> x: " << vec.x << "\t\ty: " << vec.y << "\t\tz: " << vec.z << '\n';
	}

};
