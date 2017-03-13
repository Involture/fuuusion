import math

class point :
	"""A 3D point/vector in cartesian coordinates system"""
	def __init__(self, x=0, y=0, z=0) :
		"""float * float * float -> point
		Point constructor"""
		self.x=x
		self.y=y
		self.z=z

	def __getattr__(self, attr) :
		if attr == "coords" :
			return (self.x, self.y, self.z)
		else :
			raise AttributeError

	def __str__(self) :
		"""point -> string"""
		return "{} {} {}".format(self.x,self.y,self.z)

	def __repr__(self) :
		"""point -> string"""
		return "({},{},{})".format(self.x,self.y,self.z)

	def __abs__(self) :
		"""point -> float
		Returns the euclidean norm (2-norm) of the vector"""
		return math.sqrt((self.x)**2 + (self.y)**2 + (self.z)**2)

	def __rmul__(self, scalar=1) :
		"""float * point -> point
		Returns the vector multiplied by a scalar"""
		return point(self.x * scalar, self.y * scalar, self.z * scalar)


	def normalize(self) :
		"""point -> point
		Returns a normalized version of the vector"""
		return (1/abs(self))*self

	def __add__(self, otherPoint) :
		"""point * point -> point
		Adds two points, coordinates to coordinates (vector style)"""
		return point(self.x + otherPoint.x, self.y + otherPoint.y, self.z + otherPoint.z)

	def __sub__(self,  otherPoint) :
		"""point * point -> point
		Substracts two points, coordinates to coordinates (vector style)"""
		return self + (-1) * otherPoint

	def __mul__(self, otherPoint) :
		"""point * point -> float
		Scalar product"""
		return sum(map((lambda x,y : x*y), self.coords, otherPoint.coords))


class line :
	"""A line, in 3D"""
	def __init__(self, p1=point(), p2=point()) :
		"""point * point -> line
		Line constructor, uses any two distinct points of the line"""
		self.p1 = p1
		self.p2 = p2

	def intersection(self, intersector) :
		"""line * plane -> point
		Returns the intersection between the line and the plane"""
		if intersector.normal * (self.p2 - self.p1) == 0 :
			raise "Cannot take the intersection of a line and a plane that are parallel"
		else :
			if intersector.normal.x * (self.p2 - self.p1).x != 0 :
				return (self.p2 - intersector.origin).x/(self.p2 - self.p1).x * (self.p1 - self.p2) + self.p2 #I have
			elif intersector.normal.y * (self.p2 - self.p1).y != 0 :
				return (self.p2 - intersector.origin).y/(self.p2 - self.p1).y * (self.p1 - self.p2) + self.p2 #No idea
			elif intersector.normal.z * (self.p2 - self.p1).z != 0 :
				return (self.p2 - intersector.origin).z/(self.p2 - self.p1).z * (self.p1 - self.p2) + self.p2 #What I'm doing
			#(But it should work)
			else :
				raise "You've met a terrible fate haven't you ?"


class triangle :
	"""Just an oriented 3-dimensional triangle
	Orientation is defined by the order of the points"""
	def __init__(self, p1=point(), p2=point(), p3=point()) :
		"""point * point * point -> triangle
		Triangle constructor"""
		self.p1 = p1
		self.p2 = p2
		self.p3 = p3

	def __getattr__(self, attr) :
		if attr = "points" :
			return [p1, p2, p3]
		else :
			raise AttributeError

	def trim(self, cutter) :
		"""triangle * plane -> triangle list
		Cuts the triangle so that only the part that is on the good side of the plane is conserved
		\"Good\" being defined by the normal vector"""
		numberToKeep = len(filter(lambda x : cutter.equation(x)>=0, self.points))
		if numberToKeep == 0 :
			return []

		elif numberToKeep == 1 :
			if cutter.equation(self.p2) >= 0 :		#v
				tmp = self.p1						#v
				self.p1 = self.p2					#v
				self.p2 = self.p3					#v
				self.p3 = tmp						#v
			elif cutter.equation(self.p3) >= 0 :	#Rotating the points to put the one we'll keep first
				tmp = self.p1						#^
				self.p1 = self.p3					#^
				self.p3 = self.p2					#^
				self.p2 = tmp						#^
			line12 = line(self.p1, self.p2)
			intersection2 = line12.intersection(cutter)
			line13 = line(self.p1, self.p3)
			intersection3 = line13.intersection(cutter)
			return [triangle(self.p1, intersection2, intersection3)] #keep a smaller triangle

		elif numberToKeep == 2 :
			if cutter.equation(self.p2) <= 0 :		#v
				tmp = self.p1						#v
				self.p1 = self.p2					#v
				self.p2 = self.p3					#v
				self.p3 = tmp						#v
			elif cutter.equation(self.p3) <= 0 :	#Same as above, but we put the one we don't keep first
				tmp = self.p1						#^
				self.p1 = self.p3					#^
				self.p3 = self.p2					#^
				self.p2 = tmp						#^
			line12 = line(self.p1, self.p2)
			intersection2 = line12.intersection(cutter)
			line13 = line(self.p1, self.p3)
			intersection3 = line13.intersection(cutter)
			return [triangle(intersection2, self.p2, self.p3), triangle(intersection2, self.p3, intersection3)]
		elif numberToKeep == 3 :
			return [self]



class object :
	"""A 3D object represented by a set of triangles"""
	def __init__(self, triangles=[]) :
		self.triangles = triangles

	def trim(self, cutter) :
		"""object * plane -> object
		Cuts the object so that only the part that is on the good side of the plane is conserved
		\"Good\" being defined by the normal vector"""
		newTriangles = []
		for t in self.triangles :
			newTriangles += t.trim(cutter)
		self.triangles = newTriangles


class plane :
	"""An oriented plane"""
	def __init__(self, origin=point(), normal=point(1)) :
		"""point * point -> plane
		Plane constructor, takes in argument the origin of the plane and its normal vector"""
		self.origin = origin
		try :
			self.normal = normal.normalize()
		except ZeroDivisionError :
			print("normal vector must be not null") 

	def equation(self, p) :
		return (p-self.origin) * self.normal


