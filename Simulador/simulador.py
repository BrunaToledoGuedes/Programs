import simpy
import math
import random
import argparse
from gzip import GzipFile
import sys
from scipy.special import erfc

# Assumptions:
# - All communication is assumed to be from STAs to AP.
# - AP is assumed to be fixed at coordinates at the center of the scenario.
# - Groups are determined by array groups, which is indexed by STA id and can
# change dinamically.
# - All times are in us to avoid issues with float point precision.
# - All times are in us to avoid issues with float point precision.

##
# List of the events found on the output log and their formats. Notice that not
# all events are logged by default. You may have to increase verbosity in order
# to get those.

# Entry of the received power matrix.
#  - 'idS' is the link's source node id.
#  - 'idD' is the link's destination node id.
#  - 'val' is the value of the received power in that link in dBm.
## PM idS -> idD @ val

# New packet generated at the application layer.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## + now _id_ pktId

# Transmission of an application layer packet is defered until the source node's group next slot.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'timeToGroup' is the time, in us, until the begining of the next slot assigned to the source node's group.
## D now _id_ pktId timeToGroup

# Transmission of an application layer packet will be attempted within the source node's group current slot.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'endOfSlot' is the time instant, in us, at which the source node's group current slot will end.
## G now _id_ pktId endOfSlot

# A new value has been randomly selected from the contention window for the backoff counter.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'backoffValue' is the randomly selected value for the backoff counter.
## Cw now _id_ pktId backoffValue

# Transmission of an application layer packet will be aborted due to the expiration of the current slot of the source node's group.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## A now _id_ pktId

# Medium is busy and the transmitter has begun to wait for it to become idle again.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Ms now _id_ pktId

# Transmitter was waiting for the medium to become idle and it has.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Mi now _id_ pktId

# Transmitter has begun to wait for the medium to remain idle for DIFS.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## MDs now _id_ pktId

# While waiting for DIFS, the medium has become busy.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## MDi now _id_ pktId

# Transmitter has successfully finished waiting for DIFS.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## MDo now _id_ pktId

# Transmitter has begun decrementing the backoff counter.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'backoffValue' is the current value of the backoff counter.
## Bs now _id_ pktId backoffValue

# Backoff count down has been interrupted. Might indicate both that the medium has become busy, or that the counter has zeroed.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Bi now _id_ pktId

# Backoff count down is now over.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Bo now _id_ pktId

# Transmission of the packet through the wireless link is over (including ack) and it was successful.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## S now _id_ pktId

# The wait for an ack has timed out. Transmission attempt has been unsuccessful.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Ato now _id_ pktId [ack]

# Transmission of the packet through the wireless link is over (including ack) and it was unsuccessful. The link layer has given up due to excessive number of failed retries. This application layer packet is definitely lost (at least that is how the transmitter sees it).
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## D now _id_ pktId

# The actual transmission of the packet (i.e., the insertion of the bits in the wireless link) has begun.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Ts now _id_ pktId

# The actual transmission of the packet (i.e., the insertion of the bits in the wireless link) is over.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## To now _id_ pktId

# The total amount of energy received by this node's wireless interface has been increased.
#  - 'id' is the id of the node that generated that packet
#  - 'oldValue' total amount of power in dBm received at the interface before this increase.
#  - 'newValue' total amount of power in dBm received at the interface after this increase.
#  - 'nCurrentTransmitters' number of currently active transmitters
## Ei now _id_ oldValue -> newValue [ nCurrentTransmitters ]

# The total amount of energy received by this node's wireless interface has been decreased.
#  - 'id' is the id of the node that generated that packet
#  - 'oldValue' total amount of power in dBm received at the interface before this decrease.
#  - 'newValue' total amount of power in dBm received at the interface after this decrease.
#  - 'nCurrentTransmitters' number of currently active transmitters
## Ed now _id_ oldValue -> newValue [ nCurrentTransmitters ]

# The actual reception of the packet (i.e., the retrieval of the bits in the wireless link) has begun.
#  - 'id' is the id of the node that is receiving that packet
#  - 'from' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Rs now _id_ _from_ pktId

# The actual reception of the packet (i.e., the retrieval of the bits in the wireless link) is over.
#  - 'id' is the id of the node that is receiving that packet
#  - 'from' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Ro now _id_ _from_ pktId

# The Packet Error Rate for the packet under reception has been calculated.
#  - 'id' is the id of the node that is receiving that packet
#  - 'from' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'per' the value of the packet error rate for that transmission.
## PER now _id_ _from_ pktId per

# The data frame reception has failed due to low SINR.
#  - 'id' is the id of the node that is receiving that packet
#  - 'from' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'SimTransmitters' maximum number of simultaneous transmissions during the reception attempt
## d now _id_ _from_ pktId SimTransmitters

# The data frame reception has been successful.
#  - 'id' is the id of the node that is receiving that packet
#  - 'from' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'SimTransmitters' maximum number of simultaneous transmissions during the reception attempt
## r now _id_ _from_ pktId SimTransmitters

# Before transmitting an ack, the received has just started to wait for SIFS
#  - 'id' is the id of the node that is receiving that packet
#  - 'from' is the id of the node that generated that packet
#  - 'pktId' is the id of the correspondent data packet (which starts at 0 and is incremented monotonically and independently for each node).
## MS now _id_ _from_ pktId [ack]

# The actual transmission of the ack (i.e., the insertion of the bits in the wireless link) has begun.
#  - 'id' is the id of the node that is receiving that packet
#  - 'from' is the id of the node that generated that packet
#  - 'pktId' is the id of the correspondent data packet (which starts at 0 and is incremented monotonically and independently for each node).
## Ts now _id_ _from_ pktId [ack]

# The actual transmission of the ack (i.e., the insertion of the bits in the wireless link) is over.
#  - 'id' is the id of the node that is receiving that packet
#  - 'from' is the id of the node that generated that packet
#  - 'pktId' is the id of the correspondent data packet (which starts at 0 and is incremented monotonically and independently for each node).
## To now _id_ _from_ pktId [ack]

# The actual reception of the ack (i.e., the retrieval of the bits in the wireless link) has begun.
#  - 'id' is the id of the node that generated that packet
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
## Rs now _id_ pktId [ack]

# The actual reception of the ack (i.e., the retrieval of the bits in the wireless link) is over.
#  - 'id' is the id of the node that generated the correspondent data packet
#  - 'pktId' is the id of the correspondent data packet (which starts at 0 and is incremented monotonically and independently for each node).
## Ro now _id_ pktId [ack]

# The Packet Error Rate for the ack under reception has been calculated.
#  - 'id' is the id of the node that generated the correspondent data packett
#  - 'pktId' is the id of the packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'per' the value of the packet error rate for that transmission.
## PER now _id_ pktId per [ack]

# The ack frame reception has failed due to low SINR.
#  - 'id' is the id of the node that generated the correspondent data packet
#  - 'pktId' is the id of the correspondent data packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'SimTransmitters' maximum number of simultaneous transmissions during the reception attempt
# d now _id_ pktId [ack] SimTransmitters

# The ack frame reception has been successful.
#  - 'id' is the id of the node that generated the correspondent data packet
#  - 'pktId' is the id of the correspondent data packet (which starts at 0 and is incremented monotonically and independently for each node).
#  - 'SimTransmitters' maximum number of simultaneous transmissions during the reception attempt
## r now _id_ pktId [ack] SimTransmitters


## MAC times
SLOT_TIME=52
SIFS=160
DIFS=SIFS + 2 * SLOT_TIME

# Duration of a symbol in us
SYMBOL_DURATION=40
# Number of bits per symbol (OFDM symbols, already removing coding bits)
BITS_PER_SYMBOL=26
# Data packet size in symbols
DATA_PACKET_SIZE=520 * 8 / BITS_PER_SYMBOL
# Time to transmit a data packet
DATA_PACKET_TIME=DATA_PACKET_SIZE * SYMBOL_DURATION
# ACK size in symbols
ACK_SIZE=39 * 8 / BITS_PER_SYMBOL
# Time to transmit an ack packet
ACK_PACKET_TIME=ACK_SIZE * SYMBOL_DURATION
# Time until ack timeout
ACK_TIMEOUT=SIFS + ACK_PACKET_TIME + SLOT_TIME
# Maximum number of retries for a packet in the MAC layer
RETRY_LIMIT=7
# Maximum size of the contention window
CW_MAX = 1023
# Initial size of the contention window
CW_MIN = 15
# Physical layer Constants, according to IEEE STD 802.11-2007, section 15.4.8.4
CS_THRESHOLD = -70.0

## Other environmental variables
BACKGROUND_NOISE = -95.0
ANTENNA_GAIN = 3.0
ANTENNA_HEIGHT = 1.0
TRANSMISSION_POWER = 15.0

# Helper functions
def dBm2mW(x):
	return math.pow(10.0, x/10.0)

def mW2dBm(P):
	return 10.0 * math.log10(P)

def sumdBmPower(a, b):
	return mW2dBm(dBm2mW(a) + dBm2mW(b))

def subtractdBmPower(a, b):
	return mW2dBm(dBm2mW(a) - dBm2mW(b))

class OutStream(object):

	def write(self, data):
		sys.stdout.write(data)

	def close(self):
		sys.stdout.flush()

class CompressedOutStream(object):

	def __init__(self):
		self.output = GzipFile(fileobj=OutStream(), mode="w")

	def write(self, data):
		self.output.write(data)

	def close(self):
		self.output.close()

def log(what, level=0):

	if level <= args.verbosity:
		outputStream.write(what + '\n')

# Class that handles the wireless medium common to all nodes.
class Medium:

	def __init__(self, numberOfNodes):

		self.nodeList = []
		self.powerMatrix = []
		for i in range(numberOfNodes):
			self.powerMatrix.append([0] * numberOfNodes)

	def addNode(self, node):

		id = node.getId()

		self.nodeList.append(node)

		for i in self.nodeList:

			#print(id, i.getPosX(), i.getPosY())
			if i != node:
				dist = math.sqrt((i.getPosX() - node.getPosX())**2 + (i.getPosY() - node.getPosY())**2)
				#print(dist)
				if dist == 0:
					loss = 0.0
				else:
					loss = -10.0 * math.log10(2.0 * ANTENNA_GAIN * math.pow(ANTENNA_HEIGHT, 4.0)) + 40.0 * math.log10(dist)
					#loss = -10.0 * math.log10(2.0 * ANTENNA_GAIN * math.pow(ANTENNA_HEIGHT, 4.0)) + 40.0 
			else:
				loss = 0.0
			#print(id,i.getId())
			self.powerMatrix[id][i.getId()] = TRANSMISSION_POWER - loss
			self.powerMatrix[i.getId()][id] = TRANSMISSION_POWER - loss

	def startNodeTransmission(self, node):

		id = node.getId()

		for i in self.nodeList:

			i.increaseReceivedEnergy(self.powerMatrix[id][i.getId()])

	def stopNodeTransmission(self, node):

		id = node.getId()

		for i in self.nodeList:

			i.decreaseReceivedEnergy(self.powerMatrix[id][i.getId()])

	def getPowerMatrix(self, source, dest):

		return self.powerMatrix[source][dest]

	def logPowerMatrix(self):

		for i in self.nodeList:
			for j in self.nodeList:

				if i == j:
					continue

				log('PM ' + str(i.getId()) + ' -> ' + str(j.getId()) + ' @ ' + str(self.powerMatrix[i.getId()][j.getId()]), 3)

	def logPER(self, outputFileName):
		"""
		f = open(outputFileName, 'w')

		for i in self.nodeList:
			for j in self.nodeList:

				if i == j:
					continue

				# Compute the SINR for the incoming packet
				SNR = self.powerMatrix[i.getId()][j.getId()] - BACKGROUND_NOISE

				# TODO: use a more complex error model. For now, we are just using a
				# simple mathematical model for the error in a BPSK modulation symbol.
				#symbolErrorProbability = math.erfc(math.sqrt(dBm2mW(SNR))) / 2
				symbolErrorProbability = erfc(math.sqrt(dBm2mW(SNR))) / 2
				receptionProbability = math.pow(1-symbolErrorProbability, DATA_PACKET_SIZE)
				f.write(str(i.getId()) + ' ' + str(j.getId()) + ' ' + str(1-receptionProbability) + '\n')

		f.close()
		"""

		f = open(outputFileName, 'w')

		for i in self.nodeList:
			for j in self.nodeList:
				if i == j: # or i.getId() == 0 or j.getId() == 0:
					continue
				else:

					# Compute the SINR for the incoming packet
					idPM = str(j.getId())+'->'+str(i.getId())
					SNR = (self.powerMatrix[i.getId()][j.getId()] - BACKGROUND_NOISE)
					#SNR = abs(self.powerMatrix[i.getId()][j.getId()] - BACKGROUND_NOISE)

					# TODO: use a more complex error model. For now, we are just using a
					# simple mathematical model for the error in a BPSK modulation symbol.
					#symbolErrorProbability = math.erfc(math.sqrt(dBm2mW(SNR))) / 2
					#print('===================================================')
					#print('id', idPM)
					#print('self.powerMatrix[i.getId()][j.getId()]', self.powerMatrix[i.getId()][j.getId()])
					#print('BACKGROUND_NOISE', BACKGROUND_NOISE)
					#print('SNR', SNR)
					#print('dBm2mW(SNR)',dBm2mW(SNR))
					#print('math.sqrt(dBm2mW(SNR))',math.sqrt(dBm2mW(SNR)))
					#print('erfc(math.sqrt(dBm2mW(SNR)))',erfc(math.sqrt(dBm2mW(SNR))))
					#print('erfc(math.sqrt(dBm2mW(SNR))) / 2',erfc(math.sqrt(dBm2mW(SNR))) / 2)
					symbolErrorProbability = erfc(math.sqrt(dBm2mW(SNR))) / 2
					#symbolErrorProbability = erfc(math.sqrt(SNR))/2
					receptionProbability = math.pow(1-symbolErrorProbability, DATA_PACKET_SIZE)
					#f.write(idPM + ' ' + str(j.getId()) + ' ' + str(i.getId()) + ' ' + str(receptionProbability) + '\n')   #original
					f.write(idPM + ' ' + str(j.getId()) + ' ' + str(i.getId()) + ' ' + str(self.powerMatrix[i.getId()][j.getId()]) + '\n')   #original
					#print("symbolErrorProbability",symbolErrorProbability)
					#print("receptionProbability",receptionProbability)

		f.close()

	def propagationModelFile(self, outputFileName):

		f = open(outputFileName, 'w')
		f.write(str(args.numberOfSTAs) + '\n')

		self.nodeListPM = []
		
		for j in self.nodeList:
			for i in self.nodeList:
				if i == j or i == 0 or j == 0:
					continue
				# distancia euclideana	
				dist = math.sqrt((i.getPosX() - j.getPosX())**2 + (i.getPosY() - j.getPosY())**2)
				#print("distancias -----------")
				#print(i.getId(), j.getId(), i.getPosX(),j.getPosX(),i.getPosY(),j.getPosY(),((i.getPosX() - j.getPosX())**2) + ((i.getPosY() - j.getPosY())**2))
				
				loss = -10.0 * math.log10(2.0 * ANTENNA_GAIN * math.pow(ANTENNA_HEIGHT, 4.0)) + 40.0 * math.log10(dist)
				#propagationModelFile.write(str(j.getId()) + ' ' + str(i.getId()) + ' ' + str(dist) + ' ' + str(loss) +'\n' )
				idPM = str(j.getId())+'->'+str(i.getId())
				#if (dist == 0):
				    #dist = (1/float(dist))
					#dist = math.sqrt(args.scenarioWidth**2 + args.scenarioHeight**2)
				nodePM = NodePM(idPM, j.getId(), i.getId(), dist, loss)
				self.nodeListPM.append(nodePM)
		#print("self.nodeListPM")			
		#for k in self.nodeListPM:
		#	print((k.getIdA()), str(k.getIdB()),str(k.getDist()))
		self.nodeListPMorderLoss = sorted(self.nodeListPM, key = NodePM.getId)
		#print("self.nodeListPMorderLoss")
		#for k in self.nodeListPMorderLoss:
		#	print((k.getIdA()), str(k.getIdB()),str(k.getDist()))
		for k in self.nodeListPMorderLoss:
			if k.getIdA() == 0 or k.getIdB() == 0:
				continue
			else:	
				f.write(str(k.getId()) + ' ' + str(k.getIdA()) + ' ' + str(k.getIdB()) + ' ' + str(k.getDist()) + ' ' + str(k.getLoss()) +'\n' )				
			#print(k.getId())
		f.close()

class NodePM:

	def __init__(self, id, idA, idB, dist, loss):

		self.id = id
		self.idA = idA
		self.idB = idB
		self.dist = dist
		self.loss = loss


	def getId(self):

		return(self.id)
	
	def getIdA(self):

		return(self.idA)
		
	def getIdB(self):

		return(self.idB)
		
	def getDist(self):

		return(self.dist)
		
	def getLoss(self):

		return(self.loss)



# Class that defines a node: either station or AP.
class Node:

	# Constants representing the states of a sta
	STATE_IDLE = 0
	STATE_CCA = 1
	STATE_DIFS = 2
	STATE_BACKOFF = 3
	STATE_TX = 4

	def __init__(self, env, id, posX, posY, medium, groups, ap, rate):

		self.env = env
		self.id = id
		self.posX = posX
		self.posY = posY
		self.medium = medium
		self.groups = groups;
		self.rate = rate

		self.state = self.STATE_IDLE
		self.DIFSCounter = 0
		self.backoffCounter = 0

		self.receivedEnergy = []
		self.receivedEnergy.append({'when': env.now, 'level': BACKGROUND_NOISE, 'howMany': 0})
		self.ap = ap

	def start(self):
		env.process(self.run())

	def run(self):

		# Initialize some internal state variables
		lastSuccessfullAttempt = -1
		currentPacket = -1

		# Each iteration of the next loop corresponds to a packet
		# transmission (perhaps, multiple attempts at the link layer).
		while True:

			# Check how long until the next packet. We assume that traffic follows
			# a possion distribution with the rate specified in the constructor.
			# So here we draw the interval until the next packet from an exponential
			# distribution.
			intervalToNextPacket = random.expovariate(self.rate)
			if intervalToNextPacket > 0:
				yield self.env.timeout(intervalToNextPacket)

			# Now we have a new packet to transmit.
			currentPacket = currentPacket + 1
			self.log("+", str(currentPacket))
			#print(str(self.env.now) + ' STA ' + str(self.id) + ' wants to transmit another packet.')

			# Check if we are currently at our groups slot.
			currentCycle = math.floor(env.now / (args.numberOfGroups * args.slotSize))
			currentGroup = math.floor((env.now - currentCycle * args.numberOfGroups * args.slotSize) / args.slotSize)
			#print("#######grupo#################")
			#print(self.groups[self.id])
			#print("#######node#################")
			#print(self.id)
			if currentGroup != self.groups[self.id]:

				# Not my group. Compute wait time until the next
				# slot in my group.
				if currentGroup < self.groups[self.id]:

					# Group slot is still within this cycle
					timeUntilMyGroup = (currentCycle * args.numberOfGroups * args.slotSize + self.groups[self.id] * args.slotSize) - env.now
				else:
					# Next slot is in the next cycle
					timeUntilMyGroup = ((currentCycle + 1) * args.numberOfGroups * args.slotSize + self.groups[self.id] * args.slotSize) - env.now

				self.log("D", str(currentPacket) + ' ' + str(timeUntilMyGroup))
				#print(str(self.env.now) + ' STA ' + str(self.id) + ': still not my group. Waiting ' + str(timeUntilMyGroup) + 'us until next opportunity.')

				yield env.timeout(timeUntilMyGroup)

				endOfSlot = env.now + args.slotSize
			else:

				# We are already within a slot of our group. Compute
				# when the slot is going to end.
				endOfSlot = currentCycle * args.numberOfGroups * args.slotSize + (self.groups[self.id] + 1) * args.slotSize

				self.log("G", str(currentPacket) + ' ' + str(endOfSlot), 1)
				#print(str(self.env.now) + ' STA ' + str(self.id) + ': we are at my groups slot, until ' + str(endOfSlot) + '.')

			# At this point, we are currently within the slot of our group.
			# Attempt medium access: CSMA/CA

			# Set the initial contention window size
			cw = CW_MIN

			# Zero the number of attempts for the current packet.
			attempts = 0

			# If we have not just successfully transmitted a packet, we may not
			# have to perform backoff (it depends on other conditions below).
			# Test if that is the case, and set the backoffNeeded flag accordingly.
			if lastSuccessfullAttempt == self.env.now:
				needsBackoff = True
			else:
				needsBackoff = False

			# Let's proactively choose a random backoff counter (even if we may
			# not use it later).
			self.backoffCounter = random.randint(0, cw)
			self.log('Cw', str(currentPacket) + ' ' + str(cw), 1)

			while True:

				# The slot for our group may have ended from the last point we
				# verified (e.g., last transmission attempt) to now. Check it
				# again.
				if self.env.now > endOfSlot:
					self.log("A", str(currentPacket))
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': Transmission aborted due to the end of group slot.')
					break

				# Is the medium is busy? In that case, we need
				# to wait an unspecified amount of time for the medium to be
				# idle again. Only then we can start the DIFS count down procedure.
				# As other STAs finish their transmissions, they will update our
				# receivedEnergy level and trigger the channelIdle event created
				# below.
				self.state = self.STATE_CCA
				#self.log("receivedEnergy", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))
				#if(self.receivedEnergy[-1]['level'] == -95):
					#self.log("receivedEnergy95", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))

				#print(str(self.env.now) + ' uiSTA ' + str(self.id) + ':  idle ?...' + str(self.receivedEnergy[-1]['level']) + ' > ' + str(CS_THRESHOLD))
				if self.receivedEnergy[-1]['level'] > CS_THRESHOLD:
					#print(str(self.env.now) + ' oiSTA ' + str(self.id) + ': medium became idle...' + str(self.receivedEnergy[-1]['level']))
					needsBackoff = True
					#self.log("receivedEnergyMS", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))
					self.log("Ms", str(currentPacket))
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': waiting for medium to become idle...')
					#self.log("receivedEnergychannel1", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))
					self.channelIdle = self.env.event()
					#self.log("receivedEnergychannel2", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level'])+' '+str(self.channelIdle))
					yield self.channelIdle
					#self.log("receivedEnergychannel3", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level'])+' '+str(self.channelIdle))
					#self.log("receivedEnergyMI", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))
					self.log("Mi", str(currentPacket))
					#if((self.receivedEnergy[-1]['level'])<-95):
					#	import os
					#	os.system("pause")

				# We can only proceed (backoff or transmission) if the medium has
				# been free for at least DIFS, so we wait for that to happen.
				self.state = self.STATE_DIFS
				lastDifsAttempt = env.now
				#self.log("receivedEnergyMDS", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))
				self.log("MDs", str(currentPacket))
				#print(str(self.env.now) + ' STA ' + str(self.id) + ': starting difs countdown...')
				self.difsAction = self.env.event()
				yield self.difsAction | self.env.timeout(DIFS)
				if self.env.now - lastDifsAttempt < DIFS:
					self.log("MDi", str(currentPacket))
					#print(str(self.env.now) + ' STA ' + str(self.id) + 'Medium not free for enough time (DIFS)...')
					self.state = self.STATE_IDLE
					needsBackoff = True

					continue

				self.log("MDo", str(currentPacket))
				#print(str(self.env.now) + ' STA ' + str(self.id) + ': DIFS countdown is over...')

				# Do we need to perform a backoff?
				#self.log("receivedEnergyBS", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))
				if needsBackoff == True:

					# If we got this far, then the medium has been idle for the
					# required amount of time. Now we can decrease the backoff
					# counter while it remains idle.
					self.state = self.STATE_BACKOFF
					lastBackoffAttempt = env.now
					self.log("Bs", str(currentPacket) + ' ' + str(self.backoffCounter))
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': continuing backoff countdown...')

					self.backoffAction = self.env.event()
					yield self.backoffAction | self.env.timeout(self.backoffCounter * SLOT_TIME)

					self.log("Bi", str(currentPacket))
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': Backoff count down interrupted (or done)...')

					# Either the backoff count down is over, or it was interrupted
					# because the medium became busy. Update the backoff counter
					# based on the current time in order to find what which.
					if (env.now - lastBackoffAttempt < self.backoffCounter * SLOT_TIME):
						self.backoffCounter = self.backoffCounter - math.floor((env.now - lastBackoffAttempt) / SLOT_TIME)
						self.state = self.STATE_IDLE
						continue

					self.log("Bo", str(currentPacket))
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': Backoff is over...')

				# At this point, the we probably can proceed to the transmission
				# itself. However, because of all the time we had to spend
				# before, our group slot may be over (or close). Test if the
				# transmission fits the current group slot.
				if self.env.now + DATA_PACKET_TIME > endOfSlot:

					# No, it doesn't.
					self.log("A", str(currentPacket))
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': Transmission aborted due to the end of group slot.')

					# Check if we are currently at our groups slot.
					currentCycle = math.floor(env.now / (args.numberOfGroups * args.slotSize))
					currentGroup = math.floor((env.now - currentCycle * args.numberOfGroups * args.slotSize) / args.slotSize)

					if currentGroup < self.groups[self.id]:

						# Group slot is still within this cycle
						timeUntilMyGroup = (currentCycle * args.numberOfGroups * args.slotSize + self.groups[self.id] * args.slotSize) - env.now
					else:
						# Next slot is in the next cycle
						timeUntilMyGroup = ((currentCycle + 1) * args.numberOfGroups * args.slotSize + self.groups[self.id] * args.slotSize) - env.now

					yield self.env.timeout(timeUntilMyGroup)

					break

				# Yes, it does. Proceeed to transmission.
				self.state = self.STATE_TX
				yield self.env.process(self.transmit(currentPacket))

				# Wait for ack.
				self.ackAction = self.env.event()
				yield self.ackAction | self.env.timeout(ACK_TIMEOUT)
				if self.ackAction.triggered == True:

					# Success.
					self.log("S", str(currentPacket))
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': Packet transmission completed successfully.')
					lastSuccessfullAttempt = self.env.now

					break

				else:

					# Something went wrong.
					self.log("Ato", str(currentPacket) + ' [ack]')
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': Transmission attempt failed.')

					# Check if the maximum retry limit was reached.
					attempts = attempts + 1
					if attempts > RETRY_LIMIT:
						self.log("D", str(currentPacket))
						#print(str(self.env.now) + ' STA ' + str(self.id) + ': Packet transmission failed due to retry limit.')
						break

					# Update contantion window size.
					if cw < CW_MAX:
						cw = 2 * (cw + 1) - 1

					# Choose new random backoff counter for the next attempt.
					self.backoffCounter = random.randint(0, cw)
					needsBackoff = True

					self.log('Cw', str(currentPacket) + ' ' + str(cw), 1)
					#print(str(self.env.now) + ' STA ' + str(self.id) + ': Retrying with congestion window = ' + str(cw) + '.')


			# We are finally over with the CSMA/CA proceedure. Medium becomes
			# idle.
			self.state = self.STATE_IDLE

	def transmit(self, currentPacket):

		#self.log("receivedEnergyTRANSMIT", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))
		self.log("Ts", str(currentPacket))
		#print(str(self.env.now) + ': STA ' + str(self.id) + ': Starting transmission...' +' energy: '+str(self.receivedEnergy[-1]['level']))

		# First of all, we have to update the current energy level perceived by
		# each of the other nodes. We delegate that to the medium class.
		self.medium.startNodeTransmission(self)

		#print(str(self.env.now) + ': STA ' + str(self.id) + ': Just before timeout...')

		# After that, we have to inform the AP object about the transmission
		# attempt.
		self.env.process(self.ap.receiveData(self, currentPacket))

		# Compute the duration of the transmission and wait for it to be complete
		yield self.env.timeout(DATA_PACKET_TIME)

		self.log("To", str(currentPacket))
		#print(str(self.env.now) + ' STA ' + str(self.id) + ': Ending transmission...')
		# Now the transmission is over, we should update the medium.
		self.medium.stopNodeTransmission(self)
		#self.log("receivedEnergyTRANSMITend", ' '+str(currentPacket)+' ' + str(self.id) +' '+str(self.receivedEnergy[-1]['level']))

	def getId(self):

		return(self.id)

	def getPosX(self):

		return(self.posX)

	def getPosY(self):

		return(self.posY)

	def increaseReceivedEnergy(self, increment):

		currentLevel = self.receivedEnergy[-1]['level']
		howMany = self.receivedEnergy[-1]['howMany']
		newLevel = sumdBmPower(currentLevel, increment)
		self.receivedEnergy.append({'when': self.env.now, 'level': newLevel, 'howMany': howMany + 1})

		self.log("Ei", str(currentLevel) + ' -> ' + str(newLevel) + ' [ ' + str(howMany + 1) + ']', 2)
		#print(str(self.env.now) + ' STA ' + str(self.id) + ': Increasing energy level from ' + str(currentLevel) + ' to ' + str(newLevel) + '...')

		if newLevel > CS_THRESHOLD:
			if self.state == self.STATE_DIFS and self.difsAction.triggered == False:
				self.difsAction.succeed()
			elif self.state == self.STATE_BACKOFF and self.backoffAction.triggered == False:
				self.backoffAction.succeed()

		self.cleanReceivedEnergyHistory()

	def decreaseReceivedEnergy(self, decrement):

		currentLevel = self.receivedEnergy[-1]['level']
		howMany = self.receivedEnergy[-1]['howMany']

		if howMany == 1:
			# Mitigate float point approximation errors: if we are 'removing' the
			# energy corresponding to the last still active transmitter, than,
			# instead of subtracting the power, we simply assign the noise floor.
			newLevel = BACKGROUND_NOISE
		else:
			newLevel = subtractdBmPower(currentLevel, decrement)

		self.receivedEnergy.append({'when': self.env.now, 'level': newLevel, 'howMany': howMany - 1})

		self.log("Ed", str(currentLevel) + ' -> ' + str(newLevel) + ' [ ' + str(howMany - 1) + ']', 2)
		#print(str(self.env.now) + ' STA ' + str(self.id) + ': Decreasing energy level from ' + str(currentLevel) + ' to ' + str(newLevel) + '...')

		if self.state == self.STATE_CCA:
			if newLevel <= CS_THRESHOLD and self.channelIdle.triggered == False:
				self.channelIdle.succeed()

		self.cleanReceivedEnergyHistory()

	def cleanReceivedEnergyHistory(self):

		while len(self.receivedEnergy) > 1:

			if self.receivedEnergy[1]['when'] < self.env.now - DATA_PACKET_TIME:
				self.receivedEnergy.pop(0)
			else:
				break

	def receiveData(self, source, currentPacket):

		# TODO: include a propagation delay here.

		self.log("Rs", '_' + str(source.getId()) + '_ ' + str(currentPacket))
		#print(str(self.env.now) + ' AP ' + str(self.id) + ': Starting data packet reception...')

		transmissionStart = self.env.now

		# Wait for the transmission to be concluded.
		yield self.env.timeout(DATA_PACKET_TIME)

		transmissionEnd = self.env.now

		self.log("Ro", '_' + str(source.getId()) + '_ ' + str(currentPacket))
		#print(str(self.env.now) + ' AP ' + str(self.id) + ': Ending packet reception...')

		# SNIR-based error model: use the received power with respect to the source
		# and compute how likely it is that the packet was received without errors.
		receivingPower = self.medium.getPowerMatrix(source.id, self.id)
		receptionProbability = 1
		currentStateEnd = transmissionEnd
		maxSimTransmissions = 0
		for i in range(len(self.receivedEnergy) - 1, -1, -1):

			if self.receivedEnergy[i]['when'] >= transmissionEnd:
				continue

			# Check for how long the current receivedEnergy information applies
			# to the incoming packet.
			if self.receivedEnergy[i]['when'] <= transmissionStart:
				currentStateHowLong = currentStateEnd - transmissionStart
				currentStateEnd = transmissionStart
			else:
				currentStateHowLong = currentStateEnd - self.receivedEnergy[i]['when']
				currentStateEnd = self.receivedEnergy[i]['when']

			# Compute the number of symbols affected by the current energy state.
			currentStateSymbols = currentStateHowLong / SYMBOL_DURATION
			# Compute the SINR for the incoming packet
			currentStateSINR = receivingPower - subtractdBmPower(self.receivedEnergy[i]['level'], receivingPower)

			if maxSimTransmissions < self.receivedEnergy[i]['howMany']:
				#self.log('receivedEnergyMany', '%' + str(self.receivedEnergy[i]['howMany']) + '% ' + str(maxSimTransmissions) + ' ' + str(maxSimTransmissions))
				maxSimTransmissions = self.receivedEnergy[i]['howMany']
            # TODO: use a more complex error model. For now, we are just using a
			# simple mathematical model for the error in a BPSK modulation symbol.
			#symbolErrorProbability = math.erfc(math.sqrt(dBm2mW(currentStateSINR))) / 2

			symbolErrorProbability = erfc(math.sqrt(dBm2mW(currentStateSINR))) / 2
			receptionProbability = receptionProbability * math.pow(1-symbolErrorProbability, currentStateSymbols)

			currentStateEnd = self.receivedEnergy[i]['when']

			if self.receivedEnergy[i]['when'] <= transmissionStart:
				break

		self.log('PER', '_' + str(source.getId()) + '_ ' + str(currentPacket) + ' ' + str(receptionProbability), 2)
		#print(str(self.env.now) + ' AP ' + str(self.id) + ': Estimated PER = ' + str(1-receptionProbability))
		#print(str(self.env.now) + ' node ' + str(source.id) + ': Estimated PER = ' + str(receptionProbability) + 'symbolErrorProbability' + str(symbolErrorProbability))

		# Check if packet is actually going to be received and, if so, send an ack.
					#idPM = str(j.getId())+'->'+str(i.getId())
					#SNR = (self.powerMatrix[i.getId()][j.getId()] - BACKGROUND_NOISE)
					#symbolErrorProbability = erfc(math.sqrt(dBm2mW(SNR))) / 2
					#receptionProbability = math.pow(1-symbolErrorProbability, DATA_PACKET_SIZE)
		aaa = random.random()
		#print(str(source.getId()))
		#print('if aaa > receptionProbability', 'aaa', aaa, 'receptionProbability', receptionProbability)
		#import os
		#os.system("pause")
		
		if aaa > receptionProbability:
		#if random.random() > receptionProbability:
			self.log('d', '_' + str(source.getId()) + '_ ' + str(currentPacket) + ' ' + str(maxSimTransmissions))
			#print(str(self.env.now) + ' AP ' + str(self.id) + ': Packet lost due to SINR...')
		else:
			self.log('r', '_' + str(source.getId()) + '_ ' + str(currentPacket) + ' ' + str(maxSimTransmissions))
			#print(str(self.env.now) + ' AP ' + str(self.id) + ': Packet successfully received, sending ack...')

			# Send ack.
			self.log('MS', '_' + str(source.getId()) + '_ ' + str(currentPacket) + ' [ack]')
			yield self.env.timeout(SIFS)

			self.log('Ts', '_' + str(source.getId()) + '_ ' + str(currentPacket) + ' [ack]')

			# First of all, we have to update the current energy level perceived by
			# each of the other nodes. We delegate that to the medium class.
			self.medium.startNodeTransmission(self)

			# After that, we have to inform the AP object about the transmission
			# attempt.
			self.env.process(source.receiveAck(self, currentPacket))

			# Compute the duration of the transmission and wait for it to be complete
			yield self.env.timeout(ACK_PACKET_TIME)

			# Now the transmission is over, we should update the medium.
			self.medium.stopNodeTransmission(self)
			self.log('To', '_' + str(source.getId()) + '_ ' + str(currentPacket) + ' [ack]')


	def receiveAck(self, source, currentPacket):

		self.log("Rs", str(currentPacket) + ' [ack]')
		#print(str(self.env.now) + ' STA ' + str(self.id) + ': Starting ack packet reception... to ' + str(source.id))

		transmissionStart = self.env.now

		# Wait for the transmission to be concluded.
		yield self.env.timeout(ACK_PACKET_TIME)

		transmissionEnd = self.env.now

		self.log("Ro", str(currentPacket) + ' [ack]')
		#print(str(self.env.now) + ' STA ' + str(self.id) + ': Ending ack packet reception...')

		# SNIR-based error model: use the received power with respect to the source
		# and compute how likely it is that the packet was received without errors.
		receivingPower = self.medium.getPowerMatrix(source.id, self.id)
		receptionProbability = 1
		currentStateEnd = transmissionEnd
		maxSimTransmissions = 0
		for i in range(len(self.receivedEnergy) - 1, -1, -1):

			if self.receivedEnergy[i]['when'] >= transmissionEnd:
				continue

			# Check for how long the current receivedEnergy information applies
			# to the incoming packet.
			if self.receivedEnergy[i]['when'] <= transmissionStart:
				currentStateHowLong = currentStateEnd - transmissionStart
				currentStateEnd = transmissionStart
			else:
				currentStateHowLong = currentStateEnd - self.receivedEnergy[i]['when']
				currentStateEnd = self.receivedEnergy[i]['when']

			# Compute the number of symbols affected by the current energy state.
			currentStateSymbols = currentStateHowLong / SYMBOL_DURATION
			# Compute the SINR for the incoming packet
			currentStateSINR = receivingPower - subtractdBmPower(self.receivedEnergy[i]['level'], receivingPower)

			if maxSimTransmissions < self.receivedEnergy[i]['howMany']:
				maxSimTransmissions = self.receivedEnergy[i]['howMany']

			# TODO: use a more complex error model. For now, we are just using a
			# simple mathematical model for the error in a BPSK modulation symbol.
			#symbolErrorProbability = math.erfc(math.sqrt(dBm2mW(currentStateSINR))) / 2
			symbolErrorProbability = erfc(math.sqrt(dBm2mW(currentStateSINR))) / 2
			receptionProbability = receptionProbability * math.pow(1-symbolErrorProbability, currentStateSymbols)

			currentStateEnd = self.receivedEnergy[i]['when']

			if self.receivedEnergy[i]['when'] <= transmissionStart:
				break

		self.log('PER', str(currentPacket) + '  [ack] ' + str(receptionProbability))
		#print(str(self.env.now) + ' STA ' + str(self.id) + ': Estimated PER = ' + str(receptionProbability))

		if random.random() > receptionProbability:
			self.log('d', str(currentPacket) + ' [ack] ' + str(maxSimTransmissions))
			#print(str(self.env.now) + ' STA ' + str(self.id) + ': Ack Packet lost due to SINR...')
		else:
			self.log('r', str(currentPacket) + ' [ack] ' + str(maxSimTransmissions))
			#print(str(self.env.now) + ' STA ' + str(self.id) + ': Ack Packet successfully received')
			self.ackAction.succeed()

	def log(self, type, what, level=0):

		log(type + ' ' + str(self.env.now) + ' _' + str(self.id) + '_ ' + what, level)

### Main program

# Parse command line arguments in order to set simulation parameters
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-n", "--numberOfSTAs", help="number of STAs in the simulation", type=int, default=1)
parser.add_argument("-g", "--numberOfGroups", help="number of RAW grous in the simulation", type=int, default=1)
parser.add_argument("-G", "--groupsFromFile", type=str, help="read station grouping information from file. If this option is used, any value specified with -g will be ignored.", default=None)
parser.add_argument("-S", "--slotSize", help="length of the slot of each group in us", type=int, default=50e3)
parser.add_argument("-W", "--scenarioWidth", help="width of the area used for positioning nodes in m", type=int, default=1000)
parser.add_argument("-H", "--scenarioHeight", help="height of the area used for positioning nodes in m", type=int, default=1000)
parser.add_argument("-s", "--seed", help="seed for the pseudo-random number generator", type=int, default=random.randint(0, 99999999))
parser.add_argument("-r", "--rate", help="average packet generation rate for each node in packet/us", type=float, default=10000)
parser.add_argument("-l", "--length", help="simulation length in us", type=float, default=2e7) #3e6 #2e7us = 20 segundos)
parser.add_argument("-v", "--verbosity", type=int, help="increase output log verbosity", choices=[0, 1, 2, 3, 4], default=0)
parser.add_argument("-pP", "--printPositions", type=str, help="create file with node positions", default=None)
parser.add_argument("-pE", "--printPER", type=str, help="create file with PER values for each link (considering background noise)", default=None)
parser.add_argument("-z", "--zip", help="generate zipped output", default=False, action='store_const', const=True)
parser.add_argument("-mp", "--propagationModel", help="calculates loss of path between two stations", default=None)

args = parser.parse_args()

# Create the output stream for the simulation log. Check if the user requested
# a zipped output.
if args.zip == False:
	outputStream = OutStream()
else:
	outputStream = CompressedOutStream()

# Set seed for pseudo-random number generation.
random.seed(args.seed)

# Create simulation environment
env = simpy.Environment()

# Create the medium object
medium = Medium(args.numberOfSTAs + 1)

# Create AP.
ap = Node(env, 0, args.scenarioWidth / 2.0, args.scenarioHeight / 2.0, medium, -1, None, None)
medium.addNode(ap)

# Did the user request logging nodes' positions?
if args.printPositions != None:
	positionsFile = open(args.printPositions, 'w')
	# The AP is always at the center of the scenario.
	# positionsFile.write('0 ' + str(args.scenarioWidth / 2.0) + ' ' + str(args.scenarioHeight / 2.0) + '\n')
	positionsFile.write('ap ' + str(args.scenarioWidth / 2.0) + ' ' + str(args.scenarioHeight / 2.0) + '\n')

else:
	positionsFile = None

# Check if the user specified a file with the desired grouping of the stations.
if args.groupsFromFile != None:
	i = 0
	nodeList = []
	groups = [None]
	usedCoordinates = {}
	# Yes, parse that file and generate a dictionary with the mapping between stations
	# and groups.
	groupsFile = open(args.groupsFromFile, 'r')

	# Initialize the dictionary
	coordinatesToGroups = {}

	# In this file's format, each line contains the information regarding the grouping
	# of a single node. Thus, we split the file into records consisting of one line each.
	data = groupsFile.readlines()
	for line in data:

		# Each line is of the format
		# C[XXX, YYY]
		# where:
		# - C is the cluster number (from 0 to the number of clusters minus 1)
		# - XXX is a float-point number representing the x-coordinate of the node position.
		# - YYY is a float-point number representing the y-coordinate of the node position.

		# Replace characters '[', ']' and ',' in order to facilitate parsing the
		# information.
		# print(line)
		line = line.replace('[', ' ')
		line = line.replace(']', '')
		line = line.replace(',', '')

		info = line.split()
		group = int(info[0])
		#print(info)

		posX = int(float(info[1]))
		posY = int(float(info[2]))
		coordinatesToGroups[str(posX) + '_' + str(posY)] = group

		node = Node(env, i + 1, posX, posY, medium, groups, ap, args.rate)
		medium.addNode(node)
		nodeList.append(node)
		if positionsFile != None:
			positionsFile.write(str(i) + ' ' + str(posX) + ' ' + str(posY) + '\n')

		i=i+1

		groups.append(coordinatesToGroups[str(posX) + "_" + str(posY)])

		#print(posX,posY,i)

else:
	GroupFileTXT = open('GroupFileTXT', 'w') #linha adicionada para atender ao algoritmo HMR dia 26/06/2021
	# Iterate to create stations
	nodeList = []
	groups = [None]
	usedCoordinates = {}
	for i in range(args.numberOfSTAs):
		while True:
			# TODO: change random function used below to generate real values, instead of integers.
			posX = random.randint(0, args.scenarioWidth)
			posY = random.randint(0, args.scenarioHeight)

			if str(posX) + "_" + str(posY) in usedCoordinates:
				continue
			else:
				usedCoordinates[str(posX) + "_" + str(posY)] = 1
				break

		node = Node(env, i + 1, posX, posY, medium, groups, ap, args.rate)
		medium.addNode(node)
		nodeList.append(node)

		groups.append(i % args.numberOfGroups)

		if positionsFile != None:
			#positionsFile.write(str(i % args.numberOfGroups) + ' ' + str(i) + ' ' + str(posX) + ' ' + str(posY) + '\n')
			#print(str(i+1),str(posX),str(posY))
			positionsFile.write(str(i+1) + ' ' + str(posX) + ' ' + str(posY) + '\n')

		GroupFileTXT.write(str(i % args.numberOfGroups) + ' ' + str(i+1) + ' ' + str(posX) + ' ' + str(posY) + '\n') #linha adicionada para atender ao algoritmo HMR dia 26/06/2021



if positionsFile != None:
	positionsFile.close()

#print("self.nodeList Coordinates")			
#for k in nodeList:
#	print((k.getId()), str(k.getPosX()),str(k.getPosY()))

if args.propagationModel != None:
	medium.propagationModelFile(args.propagationModel)
	
# For debug purposes
medium.logPowerMatrix()
if args.printPER != None:
	medium.logPER(args.printPER)

# Start each nodes' process
for node in nodeList:
	node.start()

env.run(until=args.length)

outputStream.close()