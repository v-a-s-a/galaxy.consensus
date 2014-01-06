class ambiguousConcordance(Exception):
  '''
  Simple exception thrown when a a sample/site contains multiple genotypes at a given threshold.
  This class is defined to make this handling more explicit than using base exceptions.
  '''
  def __init__(self, value):
     self.value = value

  def __str__(self):
     return self.value
  
class discordant(Exception):
  '''
  simple exception thrown when a collection of genotypes is not concordant.
  '''
  def __init_(self, value):
      self.value = value
  def __str__(self):
      return self.value

