# Quantum Algorithms for Shapley Value Calculation

## Authors

<a href="https://github.com/iain-burge/iain-burge">Iain Burge</a>

<a href="https://carleton.ca/scs/people/michel-barbeau/">Michel Barbeau</a>

<a href="http://www-public.imtbs-tsp.eu/~garcia_a/web/">Joaquin Garcia-Alfaro</a>

## Resources

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/python/">Python Code</a>.

<a href="https://github.com/iain-burge/QuantumShapleyValueAlgorithm/tree/main/matlab/">Matlab Code</a>.

## Summary of Results 

### Random Voting Games & Quantum Shapley Values
 
{
 "cells": [
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Random Voting Games & Quantum Shapley Values\n",
    "This notebook generates random voting games and uses a quantum algorithm to estimate the Shapley value of each player. \n",
    "It then performs some basic data analysis on the predictions."
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Import Libraries"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 46,
   "metadata": {},
   "outputs": [],
   "source": [
    "import quantumBasicVotingGame as vg\n",
    "from quantumShapEstimation import QuantumShapleyWrapper as qsw\n",
    "\n",
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "import pickle\n",
    "from tqdm.auto import tqdm"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Define Variables"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 47,
   "metadata": {},
   "outputs": [],
   "source": [
    "numTrails = 32\n",
    "\n",
    "maxEll    = 7\n",
    "\n",
    "#Defining the different conditions\n",
    "numPlayersCond    = [4,8,12]\n",
    "thresholdBitCond  = [4,5,6]\n",
    "roughVarianceCond = [1,2,2]\n"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Run Simulations"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 48,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "353ad19873484d85b191564c5d014cf5",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Trial:   0%|          | 0/32 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "0760b9d057a24e8390f8f8e0caebcd01",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "9befaa47222f48fb840ef40fa8033733",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "19bed865d4394f8eb24a31bc1013b1ce",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "36fcf2f9fa854bb3a8a4fc56be39395e",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "3f5317b05d10419e9a9952b7ec1a1e3e",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "67a529282747423b87fc6b88f3c8ddd4",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "16df4f5d03cd4c9d8c3f8af914640195",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "9949a71471964a91976aedf211b0d441",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "a8d169cf001e42c9ba3947d7aa71148c",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "94cbe7e36ef94467b673ac872e89ad01",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "d3e97c6ef17241ba84065f6858a6f4e4",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "78290f474dbe465c9f8dc32b61ac85e3",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "04af1e37e3ef445d9e9f166383df42a4",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "0c2d5faea455416b982ac2b1f4951f65",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "48ab54beca4844618ace3428c8ee2ca0",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "78f3e89d0f24451db9fe58153d675dc0",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "5f76714145ab466b8661a0d9393a4cd2",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "be97a6aef43249ba85d73247b7d9c558",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "3f9a9906dd66427080b36a2a8a508772",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "d18776efc9a343d2bff47659ba7d2456",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "96fb13f06f9044bea80d7b420aa8e755",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "34cd2f38fd034a39a8d56acddf586f94",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "63261b2ee3b9401ba4bfe4fadafa9872",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "119ab0e0a2fe4a6fa687e05c8d66fc00",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "5cef775c249c429bb06ed2738a095f51",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "bee397c2997a4308ac0bfe2f122f3462",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "b84c0dc39f42451eb44fbdc324b2ac9f",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "69849a0629b84278a3736faa0f110af6",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "6576855ac23f47919a0f3e2c8d026f70",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "ace449c98d24441cb5e72a0e6fb347ce",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "e6b3b988b5db461886c13c07a2578d92",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    },
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "e0e15151b8d546af9218dbfafa546e59",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Current Ell:   0%|          | 0/5 [00:00<?, ?it/s]"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "simulations = {}\n",
    "\n",
    "for trialNum in tqdm(range(numTrails), desc=\"Current Trial\"):\n",
    "    for ell in tqdm(range(1,maxEll), desc=\"Current Ell\"):\n",
    "        for n, thresholdBits, roughVariance in zip(\n",
    "            numPlayersCond, thresholdBitCond, roughVarianceCond\n",
    "        ):\n",
    "            trial = (n,ell,trialNum)\n",
    "\n",
    "\n",
    "            #New random game\n",
    "            threshold = 2**(thresholdBits-1)\n",
    "\n",
    "            playerVals = vg.randomVotingGame(\n",
    "                numPlayers=n,\n",
    "                thresholdBits=thresholdBits,\n",
    "                roughVariance=roughVariance\n",
    "            )\n",
    "\n",
    "            #quantum Shapley\n",
    "            qshaps = vg.quantumVotingShap(\n",
    "                threshold=threshold,\n",
    "                playerVals=playerVals,\n",
    "                ell=ell\n",
    "            )\n",
    "\n",
    "            #classical Shapley\n",
    "            cshaps = vg.classicalVotingShap(\n",
    "                threshold=threshold,\n",
    "                playerVals=playerVals,\n",
    "            )\n",
    "\n",
    "            #Store outcome\n",
    "            simulations[trial] = (qshaps, cshaps)\n",
    "\n"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Save Results"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 49,
   "metadata": {},
   "outputs": [],
   "source": [
    "with open('shapleyVoteResults.pkl', 'wb') as f:\n",
    "    pickle.dump(simulations, f)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 50,
   "metadata": {},
   "outputs": [],
   "source": [
    "# with open('shapleyVoteResults_Apr30_11-55PM.pkl', 'rb') as fp:\n",
    "#     simulations = pickle.load(fp)"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Analyze Trials "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 51,
   "metadata": {},
   "outputs": [],
   "source": [
    "def meanAbsError(qshaps, cshaps):\n",
    "    err = 0\n",
    "    for qshap, cshap in zip(qshaps, cshaps):\n",
    "        err += abs(qshap-cshap)\n",
    "    return err"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 52,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAA1gAAAFgCAYAAACmKdhBAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjQuMywgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy/MnkTPAAAACXBIWXMAAAsTAAALEwEAmpwYAAA5iElEQVR4nO3debhmVXnn/e9PUHBgpiSEwcII2mg7kBJItBMDAYEYyzZKMEYR6a6YYMQhrZjkDelE+9XuRMR26ooQCmNE4shriFhBnNICFogD4FBBkSJglTIqahju94+9Sh7LOqeeqrP3c6bv57rOdfZeez37uU+pt/vee+21UlVIkiRJkmbuAbMdgCRJkiQtFBZYkiRJktQTCyxJkiRJ6okFliRJkiT1xAJLkiRJknpigSVJkiRJPbHA0qxI8qIkn53tOCQtDkmWJqkk2892LJIWlyTnJHndbMehybHA0jZJcmCSHyX5u2n6/HmSu5N8P8ltSf5vkl+aZJyS5rdWGF2Y5NYkNyd561RFUpKnJbmv5Zw7k3wtyUmTjlnS/JDkpUnWJPlxknM2OXZ4ktVJbkmyIck/JNl7mnN9sl0XfT/Jd5N8cLr+WtgssLSt3gZ8fox+76uqhwFLgM8CH0ySQSObgneupXnp7cB6YG/gicCvAn8wTf9/azlnZ+A1wN8kOXjoIDcnHf9/Vpq7/g14HXD2Zo7tBqwElgKPAO4E/nYL53tpyz8HAbsCZ/QV6Nbymmd2mfi11ZKcANwGXDzuZ6rqbmAV8HPAHps555lJbkhyR5Irkvyn1v5zSe5KssdI30Pa3aQHtv0XJ7m23eG+KMkjRvpWklOSfAP4RrvgOSPJ+vZdX07yuG38p5A0vAOA86vqR1V1M/Ax4LFb+lB1PgzcCvxMgZXkpJY37kxyXZLfGzn2lSS/ObL/wHZH+klt//D2RP62JF9M8rSRvp9M8vok/wLcBTyyDYm+rn3XN5M8f1v/MST1p6o+2PLE9zZz7J+q6h+q6o6qugt4K/CUMc97C/AB4GeuL5LsluSj7Trm1ra9bzv23CRXbNL/lUk+0rZ3SPJXSb6d5DtJ3pnkwe3Y05KsS/KaJDcDf5tkz3b+29qTuM9402cy/EfWVkmyM/AXwCu38nM7AC8Cbqiq726my+fp7k7vDvw98A9JdmwXVJ8Ejh/p+wLgvKq6O8ly4I+BZ9M9JfsM8N5Nzv0s4DC6i6yjgV+hu7u0SzvvzyRWSXPGm4ETkjwkyT7AsXRF1rSSPCDJf6a7i/zlzXRZDzyD7knXScAZSQ5px84Ffnek73HATVX1hRbDP9Ld9d4d+CPgA0mWjPR/AbAC2AnYALwFOLaqdgJ+Gbhqy3+2pDnmV4Crx+mYZE/gt4AvbObwA+iehD0C2B/4IV3xBnABcECS/zDS/wV0OQngDXTXL08EHgXsA/zZSN+fo8tLj6DLQa8C1tFdH+1Fd71U4/wNmhkLLG2tvwTOqqp1Y/Y/PsltwA3ALwL/eXOdqurvqup7VXVPVf01sAPw6HZ4Fe1iJ8l2wPOAd7djLwH+36q6tqruAf4H8MTRp1jt+C1V9UPgbrqLnscAaZ+7acy/RdLkfZruidUddBcKa4APT9P/51vO+S5wOvCCqvrapp2q6h+r6l/bk65PAR8H/lM7/HfAce2GEnQXOBtzzu8CF1bVhVV1X1WtbjEdN3L6c6rq6paT7gHuAx6X5MFVdVNVjXWRJmluSPJ4ukLmv22h61ta/vkicBObuRndrnU+UFV3VdWdwOvphj5TVT8G3sf91zyPpRui+NH2esUK4BXtmuZOumueE0ZOfx9welX9eOSaZ2/gEVV1d1V9pqossCbAAktjS/JE4NfZujHF51fVrlX18Ko6oqqu2FynJH/Uhuvc3pLTLsCe7fBHgIOTHAAcBdxeVZe3Y48AzmyPv28DbgFCd1dnoxs2blTVJ+juFL0NWJ9k5chFlKQ5pA1l+RjwQeChdDlhN+CN03zs31rO2b2qnlhV501x7mOTXNqGzdxGVyDtCVBV/wb8C/BbSXale2r2nvbRRwDP3Zhz2mefSncRs9FozvkB8Nt0N4NuSvKPSR6zNf8OkmZPkkcB/wScWlWf2UL3l7X8s09VPb+qNmzmfA9J8n+SXJ/kDrqbSLu2G8jQ3VT+nVZQvYDuOurHdE+hHgJcMZJ7PtbaN9pQVT8a2f9fwFrg422Y8mlb/Q+gbWKBpa3xNLo7Kd9u43v/iO4C5MqZnDTd+1avphuut1tV7QrcTlco0ZLF+XR3dEbvJEN3IfN7LaFt/HlwVf3fkT4/dbemqt5SVb9IN2TwILZ8R0rS7NidbgjNW9sd2e/RDa05bvqPTa8NWf4A8FfAXi3nXEjLOc3GJ+fPBT5XVTe29huAd2+Scx5aVW8Y+eymOeeiqjqKrgj7KvA3M4lf0mS00TD/DPxlVb17S/3H9Cq6ETqHVdXOdEMP4f5rnkuBf6d7ov473H/N81264YSPHck9u7RJNTbaNPfcWVWvqqpHAs8EXpnkyJ7+Dk3DAktbYyXwC3Rjf58IvJPuXYSnz/C8O9ENo9kAbJ/kz+jeixh1Lt07XM/kpwusdwKvbY/RSbJLkudO9UVJnpzksHQTZPwA+BHdI3VJc0x7X/ObwO8n2b49TToR+NIMT/0gumHIG4B7khxL937mqA8DhwCncv/7D9ANH/zNJE9Psl2SHdvL5ftu7ouS7JVkeZKHAj8Gvo85R5oTWl7ZEdgO2Pi/5+3bsX2AT9Dd4Hlnj1+7E12hdFuS3emGMm/qXLrRNndX1WcBquo+upszZyR5+MYYk0x5DZbkGUke1Z6G3Q7ci/lnIiywNLY2XvjmjT90Fwo/2twj8K10Ed1j7q8D19MVPTeMdqiqf6FLCldW1fUj7R+iGy50XnvU/hW64TxT2ZkuQd3avut7dI/QJc1NzwaOoSuG1tK9U/CKmZywvbvwMron47fS3SW+YJM+P6R7ynUA3RDFje03ABsn19lAl6v+G1P//+kD6N7D+De6Icy/Cvz+TOKX1Js/pSt2TqN7Yv3D1gbwX4BHAn+ebm2r7yf5fg/f+WbgwXRPpC5l85P2vJtuBsJN1xp9DV0evLRd8/wz97+vvjkHtj7fBz4HvL2qLplJ8BpPfNdN80WSTwB/X1Xvmu1YJC187Wn6QVX1u1vsLEk9aVOvrwcOqapvzHY82nouQqZ5IcmT6YbrLJ/tWCQtfG3ozsl0731K0iT9PvB5i6v5ywJLc16SVXRrWZ3ahvZI0mCS/Fe6YTzvrqpPz3I4khaRJN+im/DiWbMbiWbCIYKSJEmS1BMnuZAkSZKknizIIYJ77rlnLV26dLbDkARcccUV362qJVvuOX+Zc6S5wXwjaZKmyjkLssBaunQpa9asme0wJAFJrt9yr/nNnCPNDZPKN0nOBp4BrK+qx7W2/wX8Jt0isf8KnFRVt7Vjr6WbNOVe4GVVdVFrPwY4k24dpndtsmD1ZplvpLljqpzjEEFJkqStcw7d+myjVgOPq6rH063r+FqAJAcDJwCPbZ95e1ukejvgbXRrNx4MPK/1lTTPWWBJkiRthTa75C2btH28qu5pu5cC+7bt5cB5VfXjqvom3UKxh7aftVV1XVX9O3AeLkUiLQgWWJIkSf16MfBPbXsf4IaRY+ta21TtPyPJiiRrkqzZsGHDAOFK6pMFlqQFIcnZSdYn+com7X+Y5KtJrk7yP0faX5tkbZKvJXn6SPsxrW1tktMm+TdImv+S/AlwD/Cevs5ZVSurallVLVuyZEHP4SEtCAtykgtJi9I5wFuBczc2JPk1uiE3T6iqHyd5eGsffSfi54F/TnJQ+9jbgKPo7iZ/PskFVXXNxP4KSfNWkhfRTX5xZN2/0OiNwH4j3fZtbUzTLmke8wmWpAVhc+9EAL8PvKGqftz6rG/tvhMhqVdtRsBXA8+sqrtGDl0AnJBkhyQHAAcClwOfBw5MckCSB9Hd9Llg0nFL6p8FlqSF7CDgPyW5LMmnkjy5tftOhKRtluS9wOeARydZl+RkuifoOwGrk1yV5J0AVXU1cD5wDfAx4JSqurdNiPFS4CLgWuD81lfSPOcQQUkL2fbA7sDhwJOB85M8so8TV9VKYCXAsmXLagvdJS0gVfW8zTSfNU3/1wOv30z7hcCFPYYmaQ6wwJK0kK0DPtjehbg8yX3AnvhOhCRJGohDBCUtZB8Gfg2gTWLxIOC7+E6EJEkaiE+wJC0I7Z2IpwF7JlkHnA6cDZzdpm7/d+DE9jTr6iQb34m4h/ZORDvPxncitgPO9p0ISZK0NSywJC0IU7wTAfC7U/T3nQhJktQ7CyxpHjlj9ddnOwQAXnHUQVvuJGlemyv5Bsw50mIwV3JOH/nGd7AkSZIkqScWWJIkSZLUEwssSZIkSeqJBZYkSZIk9cQCS5IkSZJ6YoElSZIkST2xwJIkSZKknlhgSZIkSVJPLLAkSZIkqScWWJIkSZLUEwssSZIkSeqJBZYkSZIk9cQCS5IkSZJ6YoElSZIkST2xwJIkSZKknlhgSZIkSVJPLLAkSZIkqScWWJIkSZLUEwssSZIkSeqJBZYkSZIk9cQCS5IkSZJ6YoElSZIkST2xwJIkSZKknlhgSZIkSVJPLLAkSZIkqScWWJIWhCRnJ1mf5CubOfaqJJVkz7afJG9JsjbJl5IcMtL3xCTfaD8nTvJvkCRJ858FlqSF4hzgmE0bk+wHHA18e6T5WODA9rMCeEfruztwOnAYcChwepLdBo1akiQtKBZYkhaEqvo0cMtmDp0BvBqokbblwLnVuRTYNcnewNOB1VV1S1XdCqxmM0WbJEnSVAYvsJJsl+QLST7a9g9IclkbmvO+JA9q7Tu0/bXt+NKRc7y2tX8tydOHjlnSwpBkOXBjVX1xk0P7ADeM7K9rbVO1S5IkjWUST7BOBa4d2X8jcEZVPQq4FTi5tZ8M3Nraz2j9SHIwcALwWLo7yW9Pst0E4pY0jyV5CPDHwJ8NdP4VSdYkWbNhw4YhvkKSJM1DgxZYSfYFfgN4V9sPcATw/tZlFfCstr287dOOH9n6LwfOq6ofV9U3gbV070ZI0nR+ATgA+GKSbwH7Alcm+TngRmC/kb77trap2n9GVa2sqmVVtWzJkiUDhC9JkuajoZ9gvZnu3Yf72v4ewG1VdU/bHx1+85OhOe347a2/Q3YkbbWq+nJVPbyqllbVUrrccUhV3QxcALywzSZ4OHB7Vd0EXAQcnWS3NrnF0a1NkiRpLIMVWEmeAayvqiuG+o5Nvs/hOtIiluS9wOeARydZl+TkabpfCFxH90T8b4A/AKiqW4C/BD7ffv6itUmSJI1l+wHP/RTgmUmOA3YEdgbOpJuta/v2lGp0+M3GoTnrkmwP7AJ8jzGH7FTVSmAlwLJly2rT45IWtqp63haOLx3ZLuCUKfqdDZzda3CSFpQkZwMbbyQ/rrXtDrwPWAp8Czi+qm5trzucCRwH3AW8qKqubJ85EfjTdtrXVdUqJM17gz3BqqrXVtW+7aLmBOATVfV84BLgOa3bicBH2vYFbZ92/BPtIugC4IQ2y+ABdOvWXD5U3JIkSVtwDj+7hMNpwMVVdSBwcdsH192TFp3ZWAfrNcArk6yle8fqrNZ+FrBHa38lLTFV1dXA+cA1wMeAU6rq3olHLUmSxJTr7o1O1rXpJF6uuyctIkMOEfyJqvok8Mm2fR2bmQWwqn4EPHeKz78eeP1wEUqSJM3IXm2yHICbgb3a9ozX3Uuygu7pF/vvv3+PIUsawmw8wZIkSVqw2isOvb0P7rIQ0vxigSVJkjRz32lD/2i/17f2Ga+7J2l+scCSJEmaudHJujadxMt196RFZCLvYEmSJC0Ubd29pwF7JllHNxvgG4Dz2xp81wPHt+4X0k3RvpZumvaToFt3L8nGdffAdfekBcMCS5IkaStMs+7ekZvp67p70iLjEEFJkiRJ6okFliRJkiT1xAJLkiRJknpigSVJkiRJPbHAkiRJkqSeWGBJkiRJUk8ssCRJkiSpJxZYkiRJktQTCyxJkiRJ6okFliRJkiT1xAJLkiRJknpigSVJkiRJPZm2wEqyXZKvTioYSYubOUfSpJhvJA1l2gKrqu4FvpZk/wnFI2kRM+dImhTzjaShbD9Gn92Aq5NcDvxgY2NVPXOwqCQtZuYcSZNivpHUu3EKrP9n8Cgk6X7mHEmTYr6R1LstFlhV9akkewFPbk2XV9X6YcOStFiZcyRNivlG0hC2OItgkuOBy4HnAscDlyV5ztCBSVqczDmSJsV8I2kI4wwR/BPgyRvv6CRZAvwz8P4hA5O0aJlzJE2K+UZS78ZZB+sBmzwu/96Yn5OkbbFNOSfJ2UnWJ/nKSNv/SvLVJF9K8qEku44ce22StUm+luTpI+3HtLa1SU7r6W+SNDd5jSOpd+MkkY8luSjJi5K8CPhH4MJhw5K0iG1rzjkHOGaTttXA46rq8cDXgdcCJDkYOAF4bPvM29uaONsBbwOOBQ4Gntf6SlqYvMaR1LtphwgmCfAWupc/n9qaV1bVh4YOTNLiM5OcU1WfTrJ0k7aPj+xeCmx8t2I5cF5V/Rj4ZpK1wKHt2Nqquq7Fc17re822/UWS5iqvcSQNZdoCq6oqyYVV9R+BD04oJkmL1MA558XA+9r2PnQF10brWhvADZu0H7a5kyVZAawA2H9/1ymV5huvcSQNZZwhglcmefKWu0lSL3rPOUn+BLgHeE9f56yqlVW1rKqWLVmypK/TSposr3Ek9W6cWQQPA56f5Hq6Vc5Dd+Pn8YNGJmmx6jXntPcqngEcWVXVmm8E9hvptm9rY5p2SQuP1ziSejfOO1grgOsnE46kxazvnJPkGODVwK9W1V0jhy4A/j7Jm4CfBw6kWwsnwIFJDqArrE4AfqePWCTNLV7jSBrKOO9gva2NT5akQc0k5yR5L/A0YM8k64DT6WYN3AFY3V1LcWlVvaSqrk5yPt3kFfcAp1TVve08LwUuArYDzq6qq3v40yTNMV7jSBrKOEMEr0zy5Kr6/ODRSNI25pyqet5mms+apv/rgddvpv1CnKZZWiy8xtGCdsbqr892CAC84qiDZjuEifIdLElzjTlH0qSYbyT1bpwC6+mDRyFJ9zPnSJoU842k3k05TXuSIwCq6nrgAVV1/cYf4BcnFaCkxcGcI2lSzDeShjTdOlh/NbL9gU2O/ekAsUha3Mw5kibFfCNpMNMVWJlie3P7kjRT5hxJk2K+kTSY6QqsmmJ7c/uSNFPmHEmTYr6RNJjpJrl4ZJIL6O7kbNym7R8weGSSFhtzjqRJMd9IGsx0Bdbyke2/2uTYpvuSNFPmHEmTYr6RNJgpC6yq+tQkA5G0uJlzJE2K+UbSkKZ7B0uSJEmStBUssCRJkiSpJ2MXWEkeMmQgkjTKnCNpUvrMN0lekeTqJF9J8t4kOyY5IMllSdYmeV+SB7W+O7T9te340r7ikDR7tlhgJfnlJNcAX237T0jy9sEjk7QomXMkTUrf+SbJPsDLgGVV9ThgO+AE4I3AGVX1KOBW4OT2kZOBW1v7Ga2fpHlunCdYZwBPB74HUFVfBH5lyKAkLWrmHEmTMkS+2R54cJLtgYcANwFHAO9vx1cBz2rby9s+7fiRSVzoWJrnxhoiWFU3bNJ07wCxSBJgzpE0OX3mm6q6kW6a92/TFVa3A1cAt1XVPa3bOmCftr0PcEP77D2t/x6bnjfJiiRrkqzZsGHDtoYnaULGKbBuSPLLQCV5YJI/Aq4dOC5Ji5c5R9Kk9JpvkuxG91TqAODngYcCx8w0yKpaWVXLqmrZkiVLZno6SQMbp8B6CXAK3V2WG4EnAn8wYEySFjdzjqRJ6Tvf/DrwzaraUFV3Ax8EngLs2oYMAuzbvov2ez+AdnwX2nBFSfPXlAsNj3h0VT1/tCHJU4B/GSYkSYucOUfSpPSdb74NHN5mJfwhcCSwBrgEeA5wHnAi8JHW/4K2/7l2/BNVVdv43ZLmiHGeYP3vMdskqQ/mHEmT0mu+qarL6CaruBL4Mt111krgNcArk6yle8fqrPaRs4A9WvsrgdO29bslzR1TPsFK8kvALwNLkrxy5NDOdNOOTivJjsCngR3a97y/qk5PcgDdHZw96F78fEFV/XuSHYBzgV+kezz+21X1rXau19JNZXov8LKqumhr/1BJc9tMc44kjWvIfFNVpwOnb9J8HXDoZvr+CHjuTL5P0twz3ROsBwEPoyuOdhr5uYPuMfaW/Bg4oqqeQDem+Zgkh7OVa0EkOZhuDYnH0r0o+vYkXmxJC89Mc44kjct8I2kwUz7BqqpPAZ9Kck5VXb+1J25jiL/fdh/YfopuLYjfae2rgD8H3kE3686ft/b3A29ta0EsB86rqh8D32yP0Q+lG68saYGYac6RpHGZbyQNaZxJLs5J8jMvXFbVEVv6YHvSdAXwKOBtwL8y5loQSTauBbEPcOnIaUc/M/pdK4AVAPvvv/8Yf5akOWqbc44kbSXzjaTejVNg/dHI9o7AbwH3TNH3p1TVvcATk+wKfAh4zNYGOK6qWkn3IinLli1zBh5p/trmnCNJW8l8I6l3WyywquqKTZr+JcnlW/MlVXVbkkuAX6KtBdGeYm1uLYh1m6wF8ZM1IprRz0haYPrIOZI0DvONpCFscZr2JLuP/OyZ5Ol0xc+WPrekPbkiyYOBo+hWR9+4FgRsfi0I+Om1IC4ATkiyQ5uB8EDA5CctUNuacyRpa5lvJA1hnCGCV9BNThG6x+bf5P6Z/6azN7CqvYf1AOD8qvpokmuA85K8DvgCP70WxLvbJBa30M0cSFVdneR84Jr2/ae0oYeSFqZtyjlJzgaeAayvqse1tt2B9wFLgW8Bx1fVrW0CnTOB44C7gBdV1ZXtMycCf9pO+7qqWtXbXyZprtnWaxxJmtI4QwQP2JYTV9WXgCdtpn2r14KoqtcDr9+WOCTNL9uac4BzgLfSrae30WnAxVX1hiSntf3XAMfSPQ0/EDiMbibTw1pBdjqwjO6i64okF1TVrdsYk6Q5bAb5RpKmNN1Cw8+e7oNV9cH+w5G0WM0051TVp5Ms3aR5OfC0tr0K+CRdgbUcOLcNQ740ya5J9m59V1fVLS2m1XTr7713a/4WSXOb1ziShjTdE6zfnOZYASYfSX0aIufsVVU3te2bgb3a9k+WhWg2Lv8wVfvPcGkIaV7zGkfSYKZbaPikSQYiaXEbOudUVW1uvZsZnM+lIaR5ymscSUMaZxbBXZK8Kcma9vPXSZxhR9Iges4532lD/2i/17f2qZZ/cFkIaRHxGkfSELZYYAFnA3cCx7efO4C/HTIoSYtanzlndPmHTZeFeGE6hwO3t6GEFwFHJ9ktyW7A0a1N0sLkNY6k3o0zTfsvVNVvjez/9yRXDRSPJG1TzknyXrpJKvZMso5uNsA3AOcnORm4nu4CCuBCuina19JN034SQFXdkuQvgc+3fn+xccILSQuS1ziSejdOgfXDJE+tqs8CJHkK8MNhw5K0iG1Tzqmq501x6MjN9C3glCnOczbdXW1JC5/XOJJ6N06B9ft0CwbvQrcQ3y3Ai4YMStKiZs6RNCnmG0m9G2eh4auAJyTZue3fMXRQkhYvc46kSTHfSBrCOLMIntoSz53Am5JcmeTo4UOTtBiZcyRNivlG0hDGmUXwxe2OztHAHsAL6F4cl6QhmHMkTYr5RlLvximw0n4fB5xbVVePtElS38w5kibFfCOpd+MUWFck+Thd8rkoyU7AfcOGJWkRM+dImhTzjaTejTOL4MnAE4HrququJHvQ1oyRpAGYcyRNivlGUu/GmUXwviRLgd9NUsBnq+pDg0cmaVEy50iaFPONpCGMM4vg24GXAF8GvgL8XpK3DR2YpMXJnCNpUsw3koYwzhDBI4D/UFUFkGQVcM2gUUlazMw5kibFfCOpd+NMcrEW2H9kfz/gG8OEI0nmHEkTY76R1Lspn2Al+f+AAnYCrk1yeds/DLh8MuFJWizMOZImxXwjaUjTDRH8q2mOVd+BSFr0zDmSJsV8I2kwUxZYVfWpzbUneSrwPODTQwUlafEx50iaFPONpCGNM8kFSZ4E/A7wXOCbwAeGDErS4mbOkTQp5htJfZvuHayD6O7iPA/4LvA+IFX1axOKTdIiYs6RNCnmG0lDmu4J1leBzwDPqKq1AEleMZGoJC1G5hxJk2K+kTSY6aZpfzZwE3BJkr9JciSQyYQlaREy50iaFPONpMFMWWBV1Yer6gTgMcAlwMuBhyd5R5KjJxSfpEXCnCNpUsw3koa0xYWGq+oHVfX3VfWbwL7AF4DXDB6ZpEXJnCNpUsw3koawxQJrVFXdWlUrq+rIoQKSpI3MOZImxXwjqS9bVWBJkiRJkqZmgSVJkiRJPbHAkiRJkqSebLHASvLsJN9IcnuSO5LcmeSOSQQnafEx50ialCHyTZJdk7w/yVeTXJvkl5LsnmR1+67VSXZrfZPkLUnWJvlSkkP6+cskzaZxnmD9T+CZVbVLVe1cVTtV1c5DByZp0TLnSJqUIfLNmcDHquoxwBOAa4HTgIur6kDg4rYPcCxwYPtZAbxjht8taQ4Yp8D6TlVdO3gkktTpPeckeUWSq5N8Jcl7k+yY5IAkl7U7x+9L8qDWd4e2v7YdX9pnLJLmlF7zTZJdgF8BzgKoqn+vqtuA5cCq1m0V8Ky2vRw4tzqXArsm2buveCTNju3H6LMmyfuADwM/3thYVR8cKihJi1qvOSfJPsDLgIOr6odJzgdOAI4Dzqiq85K8EziZ7u7xycCtVfWoJCcAbwR+eyZ/kKQ5q+9rnAOADcDfJnkCcAVwKrBXVd3U+twM7NW29wFuGPn8utZ200gbSVbQPeFi//3338bQJE3KOAXWzsBdwOjK5gVYYEkawhA5Z3vgwUnuBh5Cd/FyBPA77fgq4M/pCqzlbRvg/cBbk6SqagbfL2lu6jvfbA8cAvxhVV2W5EzuHw7YnbyqkmxVPqmqlcBKgGXLlpmLpDluiwVWVZ00iUAkCfrPOVV1Y5K/Ar4N/BD4ON1d5duq6p7WbeNdYxi5o1xV9yS5HdgD+O7oeb2jLM1/A1zjrAPWVdVlbf/9dAXWd5LsXVU3tSGA69vxG4H9Rj6/b2uTNI9tscBKsiPdkJnHAjtubK+qFw8Yl6RFqu+c02brWk43dOc24B+AY2Yap3eUpfmv73xTVTcnuSHJo6vqa8CRwDXt50TgDe33R9pHLgBemuQ84DDg9pGhhJLmqXEmuXg38HPA04FP0d1duXPIoCQtan3nnF8HvllVG6rqbrqhP0+he5l8402m0bvGP7mj3I7vAnxvBt8vae4a4hrnD4H3JPkS8ETgf9AVVkcl+QZdTnpD63shcB2wFvgb4A9m+N2S5oBx3sF6VFU9N8nyqlqV5O+BzwwdmKRFq++c823g8CQPoRsieCSwBrgEeA5wHj97R/lE4HPt+Cd8/0pasHq/xqmqq4Blmzl05Gb6FnDKTL5P0twzzhOsu9vv25I8ju5u7sOHC0nSItdrzmnvQrwfuBL4Ml3eWwm8BnhlkrV071id1T5yFrBHa38lm7ygLmlB8RpHUu/GeYK1sr3D8P/Q3dl9GPBng0YlaTHrPedU1enA6Zs0Xwccupm+PwKeO5PvkzRveI0jqXfjzCL4rrb5KeCRw4YjabEz50iaFPONpCFscYhgkr2SnJXkn9r+wUlOHj40SYuROUfSpJhvJA1hnHewzgEuAn6+7X8dePlA8UjSOZhzJE3GOZhvJPVsnAJrz6o6H7gPuoU3gXsHjUrSYmbOkTQp5htJvRunwPpBkj2AAkhyOHD7oFFJWszMOZImxXwjqXfjzCL4SrqZdX4hyb8AS+jWhpGkIZhzJE2K+UZS78aZRfDKJL8KPBoI8LWqunsLH5OkbWLOkTQp5htJQ5iywEry7CkOHZSEqvrgQDFJWoTMOZImxXwjaUjTPcF6P3BV+4Huzs5GBZh8JPXJnCNpUsw3kgYzXYH1bOAE4PHAR4D3VtXaiUQlaTEy50iaFPONpMFMOYtgVX24qk4AfhX4V+Cvk3y2jVXeoiT7JbkkyTVJrk5yamvfPcnqJN9ov3dr7UnyliRrk3wpySEj5zqx9f9GkhNn9BdLmpNmmnMkaVzmG0lDGmea9h/RTVl6B/AwYMcxz30P8KqqOhg4HDglycHAacDFVXUgcHHbBzgWOLD9rADeAV1BBpwOHAYcCpy+sSiTtCBta86RpK1lvpHUu+kmuTiC7vH5ocA/A2dW1ZpxT1xVNwE3te07k1wL7AMsB57Wuq0CPgm8prWfW1UFXJpk1yR7t76rq+qWFtdq4BjgvWP/lZLmvJnmHEkal/lG0pCmewfrn4EvAZ8FdgBemOSFGw9W1cvG/ZIkS4EnAZcBe7XiC+BmYK+2vQ9ww8jH1rW2qdo3/Y4VdE++2H///ccNTdLc0VvOkaQtMN9IGsx0BdZJfXxBkocBHwBeXlV3JPdP1FNVlaT6+J6qWgmsBFi2bFkv55Q0Ub3kHEkag/lG0mCmLLCqatVMT57kgXTF1XtG1pT4TpK9q+qmNgRwfWu/Edhv5OP7trYbuX9I4cb2T840NklzSx85R5LGYb6RNKRxJrnYJukeVZ0FXFtVbxo5dAGwcSbAE+mmR93Y/sI2m+DhwO1tKOFFwNFJdmuTWxzd2iRJkiRpTpluiOBMPQV4AfDlJFe1tj8G3gCcn+Rk4Hrg+HbsQuA4YC1wF+3xfVXdkuQvgc+3fn+xccILSZIkSZpLBiuwquqz/PTK6KOO3Ez/Ak6Z4lxnA2f3F50kSZIk9W+6adr/NzDlZBHOsCOpT+YcSZNivpE0pOmeYLkehKRJMudImhTzjaTBDDqLoCSNy5wjaVLMN5KGtMV3sJIsAV4DHAzsuLG9qo4YMC5Ji5Q5R9KkmG8kDWGcadrfA1wLHAD8d+Bb3D+jnyT1zZwjaVLMN5J6N06BtUdVnQXcXVWfqqoXA97ZkTQUc46kSTHfSOrdONO0391+35TkN4B/A3YfLiRJi5w5R9KkmG8k9W6cAut1SXYBXgX8b2Bn4BWDRiVpMes95yTZFXgX8Di6qZlfDHwNeB+wlG5Y0PFVdWuSAGfSLXx+F/CiqrpyJt8vac7yGkdS77ZYYFXVR9vm7cCvDRuOpMVuoJxzJvCxqnpOkgcBDwH+GLi4qt6Q5DTgNLqX3Y8FDmw/hwHvaL8lLTBe40gawhbfwUqyqt393bi/W5KzB41K0qLVd85pd6d/BTgLoKr+vapuA5YDG6dqXgU8q20vB86tzqXArkn23tbvlzR3eY0jaQjjTHLx+HYxAkBV3Qo8abCIJC12feecA4ANwN8m+UKSdyV5KLBXVd3U+twM7NW29wFuGPn8utb2U5KsSLImyZoNGzbMIDxJs8hrHEm9G+cdrAck2a0lHZLsPubnJGlb9J1ztgcOAf6wqi5LcibdcMCfqKpKUltz0qpaCawEWLZs2VZ9VovXGau/Ptsh/MQrjjpotkOYC7zGkdS7cZLIXwOfS/IPQIDnAK8fNCpJi1nfOWcdsK6qLmv776crsL6TZO+quqkNAVzfjt8I7Dfy+X1bm6SFx2scSb3b4hDBqjoXeDbwHeAm4NlV9e6hA5O0OPWdc6rqZuCGJI9uTUcC1wAXACe2thOBj7TtC4AXpnM4cPvIUEJJC4jXOJKGMO5j8AfS3dnZuC1JQ+o75/wh8J42g+B1wEl0N5jOT3IycD1wfOt7Id0U7Wvppmk/qYfvlzR3eY0jqVfjzCJ4KvAeYE/g4cDfJfnDoQOTtDgNkXOq6qqqWlZVj6+qZ1XVrVX1vao6sqoOrKpfr6pbWt+qqlOq6heq6j9W1ZqZ/1WS5iKvcSQNYZwnWCcDh1XVDwCSvBH4HN2CfJLUN3OOpEkx30jq3TjTtAe4d2T/Xu5/lC5JfTPnSJqUQfJNku3ashAfbfsHJLksydok72vDlUmyQ9tf244vnel3S5p94zzB+lvgsiQfavvPoi3YKUkDMOdImpSh8s2pwLXAzm3/jcAZVXVeknfSPTl7R/t9a1U9KskJrd9v9/D9kmbRtE+wkjwAuJTuJe9b2s9JVfXm4UOTtNiYcyRNylD5Jsm+wG8A72r7AY6gWyICYBVdIQewvO3Tjh/Z+kuax6Z9glVV9yV5W1U9CbhyQjFJWqTMOZImZcB882bg1cBObX8P4LaquqftrwP2adv7ADe0eO5Jcnvr/93REyZZAawA2H///XsMVdIQxnkH6+Ikv+UdFUkTYs6RNCm95pskzwDWV9UVfZxvo6pa2WZCXbZkyZI+Ty1pAOO8g/V7wCuBe5P8qLVVVe08zWckaVuZcyRNSt/55inAM5McB+xI9w7WmcCuSbZvT7H2BW5s/W8E9gPWJdke2AX43jZ+t6Q5YotPsKpqp6p6QFU9sG3v5IWOpKGYcyRNSt/5pqpeW1X7VtVS4ATgE1X1fOAS4Dmt24nAR9r2BW2fdvwTVVXb+v2S5oZxnmCR5NnAU4ECPlNVHx4yKEmLmzlH0qRMKN+8BjgvyeuAL3D/TIVnAe9OspZuko0TBvhuSRO2xQIryduBRwHvbU0vSXJUVZ0yaGSSFiVzjqRJGTLfVNUngU+27euAQzfT50fAc2f6XZLmlnGeYB0B/IeNj6yTrAKuHjQqSYuZOUfSpJhvJPVunFkE1wKjc4Lu19okaQjmHEmTYr6R1LtxnmDtBFyb5HK68cmHAmuSXABQVc8cMD5Ji485R9KkmG8k9W6cAuvPBo9Cku5nzpE0KeYbSb3bYoFVVZ+aRCCSBOYcSZNjvpE0hCkLrCSfraqnJrmT7rH5Tw7hop+SembOkTQp5htJQ5qywKqqp7bfO00uHEmLlTlH0qSYbyQNaYuzCCY5PMlOI/s7JTls2LAkLVbmHEmTYr6RNIRxpml/B/D9kf0ftDZJGoI5R9KkmG8k9W6cAisbF+ADqKr7GG/2QUnaFuYcSZNivpHUu3EKrOuSvCzJA9vPqcB1QwcmadEy50iaFPONpN6NU2C9BPhl4EZgHXAYsGLIoCQtauYcSZNivpHUu3HWwVoPnDCBWCTJnCNpYsw3koYwziyCByW5OMlX2v7jk/zp8KFJWozMOZImxXwjaQjjvMj5N8B/A/4PQFV9KcnfA68bMjBpUs5Y/fXZDgGAVxx10GyHMFeYcyRNivlGUu/GeQfrIVV1+SZt9wwRjCRhzpE0OeYbSb0bp8D6bpJfAAogyXOAmwaNStJi1nvOSbJdki8k+WjbPyDJZUnWJnlfkge19h3a/tp2fOkM/xZJc5vXOJJ6N06BdQrdo/PHJLkReDndrDuSNIQhcs6pwLUj+28EzqiqRwG3Aie39pOBW1v7Ga2fpIXLaxxJvdtigVVV11XVrwNLgMcAvwo8dejAJC1OfeecJPsCvwG8q+0HOAJ4f+uyCnhW217e9mnHj2z9JS1AXuNIGsKUBVaSnZO8NslbkxwF3AWcCKwFjp9UgJIWhwFzzpuBVwP3tf09gNuqauN7FuuAfdr2PsANAO347a3/5uJdkWRNkjUbNmyYQXiSJs1rHElDmm4WwXfTDZ35HPBfgT8BAvznqrpq+NAkLTK955wkzwDWV9UVSZ7WT5idqloJrARYtmxZ9XluSYPzGkfSYKYrsB5ZVf8RIMm76F763L+qfjSRyCQtNkPknKcAz0xyHLAjsDNwJrBrku3bU6p9gRtb/xuB/YB1SbYHdgG+N4PvlzQ3eY0jaTDTvYN198aNqroXWGfikTSg3nNOVb22qvatqqXACcAnqur5wCXAc1q3E4GPtO0L2j7t+CeqyqdT0sLjNY6kwUz3BOsJSe5o2wEe3PYDVFXtPHh0khaTSeac1wDnJXkd8AXgrNZ+FvDuJGuBW+iKMkkLj9c4kgYzZYFVVdtNMhBJi9vQOaeqPgl8sm1fBxy6mT4/Ap47ZBySZp/XOJKGNM46WJIkSZKkMQxWYCU5O8n6JF8Zads9yeok32i/d2vtSfKWJGuTfCnJISOfObH1/0aSEzf3XZIkSZI0Fwz5BOsc4JhN2k4DLq6qA4GL2z7AscCB7WcF8A7oCjLgdOAwuuE8p28syiRJkiRprhmswKqqT9O9JD5qObCqba8CnjXSfm51LqWbQnlv4OnA6qq6papuBVbzs0WbJEmSJM0Jk34Ha6+quqlt3wzs1bb3AW4Y6beutU3V/jOSrEiyJsmaDRs29Bu1JEmSJI1h1ia5aGvL9La+TFWtrKplVbVsyZIlfZ1WkiRJksY26QLrO23oH+33+tZ+I7DfSL99W9tU7ZIkSZI050y6wLoA2DgT4InAR0baX9hmEzwcuL0NJbwIODrJbm1yi6NbmyRJkiTNOVMuNDxTSd4LPA3YM8k6utkA3wCcn+Rk4Hrg+Nb9QuA4YC1wF3ASQFXdkuQvgc+3fn9RVZtOnCFJkiRJc8JgBVZVPW+KQ0dupm8Bp0xxnrOBs3sMTZIkSZIGMWuTXEiSJEnSQmOBJUmSJEk9scCSJEmSpJ5YYEmSJElSTyywJEmSepBkvySXJLkmydVJTm3tuydZneQb7fdurT1J3pJkbZIvJTlkdv8CSX2wwJIkSerHPcCrqupg4HDglCQHA6cBF1fVgcDFbR/gWODA9rMCeMfkQ5bUNwssSZKkHlTVTVV1Zdu+E7gW2AdYDqxq3VYBz2rby4Fzq3MpsGuSvScbtaS+WWBJkiT1LMlS4EnAZcBeVXVTO3QzsFfb3ge4YeRj61rbpudakWRNkjUbNmwYLmhJvbDAkiRJ6lGShwEfAF5eVXeMHquqAmprzldVK6tqWVUtW7JkSY+RShqCBZYkSVJPkjyQrrh6T1V9sDV/Z+PQv/Z7fWu/Edhv5OP7tjZJ85gFliRJUg+SBDgLuLaq3jRy6ALgxLZ9IvCRkfYXttkEDwduHxlKKGme2n62A5AkSVogngK8APhykqta2x8DbwDOT3IycD1wfDt2IXAcsBa4CzhpotFKGoQFliRJUg+q6rNApjh85Gb6F3DKoEFJmjiHCEqSJElSTyywJEmSJKknFliSJEmS1BMLLEkLWpL9klyS5JokVyc5tbXvnmR1km+037u19iR5S5K1Sb6U5JDZ/QskSdJ8YoElaaG7B3hVVR0MHA6ckuRg4DTg4qo6ELi47QMcCxzYflYA75h8yJIkab6ywJK0oFXVTVV1Zdu+E7gW2AdYDqxq3VYBz2rby4Fzq3MpsOvGBUIlSZK2xAJL0qKRZCnwJOAyYK+RBT1vBvZq2/sAN4x8bF1rkyRJ2iILLEmLQpKHAR8AXl5Vd4wea2vR1Faeb0WSNUnWbNiwocdIJUnSfOZCw5IWvCQPpCuu3lNVH2zN30myd1Xd1IYArm/tNwL7jXx839b2U6pqJbASYNmyZVtVnKlfZ6z++myH8BOvOOqg2Q5BkjTLLLAkLWhJApwFXFtVbxo5dAFwIvCG9vsjI+0vTXIecBhw+8hQQknSIjdXbup4Q2fussCStNA9BXgB8OUkV7W2P6YrrM5PcjJwPXB8O3YhcBywFrgLOGmi0UqSpHnNAkvSglZVnwUyxeEjN9O/gFMGDUqSJC1YTnIhSZIkST2xwJIkSZKknlhgSZIkSVJPLLAkSZIkqScWWJIkSZLUEwssSZIkSeqJBZYkSZIk9cQCS5IkSZJ6YoElSZIkST2xwJIkSZKknlhgSZIkSVJPLLAkSZIkqScWWJIkSZLUk+1nOwAtTGes/vpshwDAK446aLZDkCRJ0iLiEyxJkiRJ6okFliRJkiT1xAJLkiRJknpigSVJkiRJPbHAkiRJkqSeWGBJkiRJUk8ssCRJkiSpJxZYkiRJktQTCyxJkiRJ6sn2sx2AJGnuOWP112c7hJ94xVEHzXYIkiSNzQJLkiRJs2qu3NTxho764BBBSZIkSeqJT7DmEe/uSJIkSXPbvHmCleSYJF9LsjbJabMdj6SFy3wjaZLMOdLCMi8KrCTbAW8DjgUOBp6X5ODZjUrSQmS+kTRJ5hxp4ZkvQwQPBdZW1XUASc4DlgPXzPTEDruTtIkFn2/AnCPNIQs+55hvtNikqmY7hi1K8hzgmKr6L23/BcBhVfXSkT4rgBVt99HA1yYY4p7Adyf4fTNhrMMw1qk9oqqWTPD7ZmScfNPaZyvn+N+1YcynWGF+xTvJWOdVvgGvcXpmrMMw1qltNufMlydYW1RVK4GVs/HdSdZU1bLZ+O6tZazDMNbFZ7Zyznz6z89YhzOf4p1Psc5VXuOMx1iHYaxbb168gwXcCOw3sr9va5OkvplvJE2SOUdaYOZLgfV54MAkByR5EHACcMEsxyRpYTLfSJokc460wMyLIYJVdU+SlwIXAdsBZ1fV1bMc1qhZeWy/jYx1GMa6QJhvemWsw5lP8c6nWCfOnNMrYx2GsW6leTHJhSRJkiTNB/NliKAkSZIkzXkWWJIkSZLUEwusGUhydpL1Sb4y27FsSZL9klyS5JokVyc5dbZjmkqSHZNcnuSLLdb/PtsxTSfJdkm+kOSjsx3LliT5VpIvJ7kqyZrZjkfjM98MY77lG5g/Ocd8M7/Nl5wzn/INzL+cM1/yDcytnOM7WDOQ5FeA7wPnVtXjZjue6STZG9i7qq5MshNwBfCsqprxSvF9SxLgoVX1/SQPBD4LnFpVl85yaJuV5JXAMmDnqnrGbMcznSTfApZV1XxZMFCN+WYY8y3fwPzJOeab+W2+5Jz5lG9g/uWc+ZJvYG7lHJ9gzUBVfRq4ZbbjGEdV3VRVV7btO4FrgX1mN6rNq8732+4D28+cvBOQZF/gN4B3zXYsWtjMN8OYT/kGzDmanPmSc+ZTvoH5lXPMN9vOAmsRSrIUeBJw2SyHMqX2SPoqYD2wuqrmaqxvBl4N3DfLcYyrgI8nuSLJitkORguf+aZ3b2b+5BzzjSZqPuQbmFc5583Mn3wDcyjnWGAtMkkeBnwAeHlV3THb8Uylqu6tqifSrWh/aJI5NzwhyTOA9VV1xWzHshWeWlWHAMcCp7QhINIgzDf9moc5x3yjiZkv+QbmR86Zh/kG5lDOscBaRNpY3w8A76mqD852POOoqtuAS4BjZjmUzXkK8Mw25vc84Igkfze7IU2vqm5sv9cDHwIOnd2ItFCZbwYxr3KO+UaTMh/zDcz5nDOv8g3MrZxjgbVItJcqzwKurao3zXY800myJMmubfvBwFHAV2c1qM2oqtdW1b5VtRQ4AfhEVf3uLIc1pSQPbS8Ak+ShwNHAnJ4dSvOT+WYY8ynnmG80KfMp38D8yTnzKd/A3Ms5FlgzkOS9wOeARydZl+Tk2Y5pGk8BXkB3B+Kq9nPcbAc1hb2BS5J8Cfg83fjkOT896DywF/DZJF8ELgf+sao+NssxaUzmm8GYb4Zhvpnn5lHOmU/5Bsw5Q5lTOcdp2iVJkiSpJz7BkiRJkqSeWGBJkiRJUk8ssCRJkiSpJxZYkiRJktQTCyxJkiRJ6okFliRJkiT1xAJLkiRJknpigaU5qa3I/dYkh892LJIWPnOOpEkx3yx8Fliaq14C7Ag8dbYDkbQomHMkTYr5ZoGzwNJcdQzwNeCqWY5D0uJgzpE0KeabBc4CS3NOkh2B7YBDgE/NcjiSFjhzjqRJMd8sDhZYmosOpEs+X62qu2c7GEkLnjlH0qSYbxaB7Wc7AGkzlgAHActnOxBJi4I5R9KkmG8WAZ9gaS76eeADwAOS7DbbwUha8Mw5kibFfLMIWGBpTkmyPd245J8D3gncO7sRSVrIzDmSJsV8s3ikqmY7BkmSJElaEHyCJUmSJEk9scCSJEmSpJ5YYEmSJElSTyywJEmSJKknFliSJEmS1BMLLEmSJEnqiQWWJEmSJPXk/wf/55X5eqIR6wAAAABJRU5ErkJggg==",
      "text/plain": [
       "<Figure size 864x360 with 3 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.rcParams['figure.figsize'] = [12, 5]\n",
    "fig, ax = plt.subplots(1, len(numPlayersCond))\n",
    "\n",
    "#We're looking to find reciprocal mean abs error per trial\n",
    "#For each trial with n players\n",
    "for i, n in enumerate(numPlayersCond):\n",
    "    #Orient data\n",
    "    resultsX = []\n",
    "    resultsY = []\n",
    "    resultErr = []\n",
    "    for ell in range(1, maxEll):\n",
    "        trialOutcomes = []\n",
    "\n",
    "        for trialNum in range(numTrails):\n",
    "            qshaps, cshaps = simulations[(n,ell,trialNum)]\n",
    "            trialOutcomes.append(\n",
    "                meanAbsError(qshaps, cshaps)\n",
    "            )\n",
    "        \n",
    "        trialOutcomes = np.array(trialOutcomes)\n",
    "        resultsX.append(ell)\n",
    "        resultsY.append(trialOutcomes.mean())\n",
    "        resultErr.append(trialOutcomes.std())\n",
    "\n",
    "        # resultsX += len(trialOutcomes) * [ell]\n",
    "        # resultsY += trialOutcomes\n",
    "    \n",
    "    ax[i].set_title(f\"{n} Players\")#, Threshold: {2**thresholdBitCond[i]}\")\n",
    "    ax[i].bar(\n",
    "        np.array(resultsX), \n",
    "        1/np.array(resultsY),\n",
    "        # yerr=resultErr,\n",
    "        align='center',\n",
    "        alpha=0.5,\n",
    "        ecolor='black',\n",
    "        capsize=10,\n",
    "    )\n",
    "    ax[i].set_xlabel(r\"$\\ell$\")\n",
    "    ax[i].set_ylabel(r\"Reciprocal Mean Absolute Error\")\n",
    "\n",
    "plt.tight_layout()\n",
    "plt.show()"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "venvThesis",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.2"
  },
  "orig_nbformat": 4
 },
 "nbformat": 4,
 "nbformat_minor": 2
}


## References

If using this code for research purposes, please cite:

Iain Burge, Michel Barbeau and Joaquin Garcia-Alfaro. Quantum Algorithms for Shapley Value Calculation. *To appear*. May 2023.

```
@inproceedings{burge-barbeau-alfaro2023Shapley,
  title={Quantum Algorithms for Shapley Value Calculation},
  author={Burge, Iain and Barbeau, Michel and Garcia-Alfaro, Joaquin},
  booktitle={To appear},
  pages={1--8},
  year={2023},
  month={May},
}
```


 
