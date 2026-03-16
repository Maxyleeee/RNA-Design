"""
Debug: Try creating agent with verbose error output.
"""
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import warnings
warnings.filterwarnings('ignore')
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

from learna_tools.learna.environment import RnaDesignEnvironment, RnaDesignEnvironmentConfig

env_config = RnaDesignEnvironmentConfig(mutation_threshold=5, reward_exponent=9.33503385734547, state_radius=32)
env_config.use_conv = True
env_config.use_embedding = True
env = RnaDesignEnvironment([(1, '((((....))))')], env_config)

print(f"env.states: {env.states}")
print(f"env.actions: {env.actions}")
print()

# Try creating agent with full traceback
from tensorforce.agents import PPOAgent
import traceback

print("Creating agent with full tracing...")
try:
    agent = PPOAgent(
        states=env.states,
        actions=env.actions,
        network='auto',
        max_episode_timesteps=500,
        batch_size=126,
        learning_rate=0.0005991629320464973,
        likelihood_ratio_clipping=0.3,
        entropy_regularization=6.762991409135427e-05,
    )
    print("Agent created successfully!")
    model = agent.model
    print(f"Model is_initialized: {getattr(model, 'is_initialized', 'N/A')}")
    print(f"Has timestep_completed: {hasattr(agent, 'timestep_completed')}")
    print(f"Has model.act: {hasattr(model, 'act')}")
    
    # try act
    state = env.reset()
    action = agent.act(states=state)
    print(f"Agent.act SUCCESS: {action}")
except Exception as e:
    print(f"EXCEPTION: {type(e).__name__}: {str(e)}")
    traceback.print_exc()
