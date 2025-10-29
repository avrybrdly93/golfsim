let pyodideReady = loadPyodide({ indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.0/full/" });

async function runSim() {
  let pyodide = await pyodideReady;
  await pyodide.loadPackage("numpy");

  let speed = parseFloat(document.getElementById("speed").value);
  let angle = parseFloat(document.getElementById("angle").value) * Math.PI / 180;
  let spin = parseFloat(document.getElementById("spin").value);

    await pyodide.runPythonAsync(await (await fetch("sim.py")).text());
    //let pyCode = (await fetch("sim.py")).text();
    //let awaitPyCode = await pyodide.runPythonAsync(await(pyCode));
    //console.log(awaitPyCode);
    //console.log(pyCode);
    // Run the simulation with the current input values
    let code = `
    import math
    speed = ${speed}
    angle = ${angle}
    spin_rpm = ${spin}
    vel0 = (speed*math.cos(angle), 0, speed*math.sin(angle))
    omega_vec = (0, spin_rpm*2*math.pi/60, 0)
    traj = simulate_trajectory((0,0,0), vel0, omega_vec)
    traj
    `;
    let traj = await pyodide.runPythonAsync(code);
  //console.log("TRAJ " + traj);
  // Load the simulation code into Pyodide
  drawTrajectory(traj);
  //animateTrajectory(traj);
}

function drawTrajectory(traj) {
  let canvas = document.getElementById("simCanvas");
  let ctx = canvas.getContext("2d");

  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.beginPath();
  ctx.moveTo(0, canvas.height);

  drawMarkers(canvas, ctx, traj);

  console.log("Math max: ", traj, traj.map(p => p[0]), Math.max(...traj.map(p => p[0])))
  let scaleX = canvas.width; // Math.max(...traj.map(p => p[0]));
  let scaleY = canvas.height; // Math.max(...traj.map(p => p[2]));
  console.log("scales: ", scaleX, scaleY);

for (let [x, y, z] of traj) {
    //ctx.lineTo(x * scaleX, canvas.height - z * scaleY);
    scaleX = canvas.width / x; // Math.max(...traj.map(p => p[0]));
    scaleY = canvas.height / z; // Math.max(...traj.map(p => p[2]));
    ctx.lineTo(x * 100, canvas.height - z * 100);
    let xPrint = x * scaleX;
    let yPrint = canvas.height - z * scaleY;
    //console.log("This is x:" + x, "This is z: " + z);
}
  
  ctx.strokeStyle = "blue";
  ctx.stroke();
  //console.log("Trajectory:", traj);
  //console.log("Canvas Height:", canvas.height);
  //console.log("ScaleX:", scaleX, "ScaleY:", scaleY);
}

function drawMarkers(canvas, ctx, traj) {
    ctx.strokeStyle = "lightgray";
    ctx.fillStyle = "black";
    ctx.font = "12px Arial";
    ctx.textAlign = "center";

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.beginPath();
    ctx.moveTo(0, canvas.height);
    let partial = [];
    let i = 0;
    // Only plot up to current point i
    if (i < traj.length) {
        i++;
        partial += (traj.slice(0, i));
        requestAnimationFrame(drawFrame);
    }


    // Find max so far
    let maxX = Math.max(...partial.map(p => p[0]));
    let maxZ = Math.max(...partial.map(p => p[2]));


    let maxMeters = Math.floor(maxX / 100) * 100;

    for (let dist = 100; dist <= maxMeters; dist += 100) {
    let xCanvas = dist * scaleX;

    // Draw vertical line
    ctx.beginPath();
    ctx.moveTo(xCanvas, canvas.height);
    ctx.lineTo(xCanvas, 0);
    ctx.stroke();

    // Label at bottom
    ctx.fillText(`${dist} m`, xCanvas, canvas.height - 5);
}

}

function animateTrajectory(traj) {
  let canvas = document.getElementById("simCanvas");
  let ctx = canvas.getContext("2d");

  ctx.clearRect(0, 0, canvas.width, canvas.height);

  console.log(traj);//Math.max(...traj.map(p => p[0])));
  let scaleX = canvas.width / Math.max(...traj.map(p => p[0]));
  let scaleY = canvas.height / Math.max(...traj.map(p => p[1]));

  let i = 0; // index of current point

  function drawFrame() {
    // Clear screen
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Draw path so far
    ctx.beginPath();
    ctx.moveTo(0, canvas.height);
    for (let j = 0; j <= i; j++) {
      let [x, z] = traj[j];
      //console.log("Traj[j] is: " + traj[j]);
      ctx.lineTo(x * scaleX, canvas.height - z * scaleY);
    }
    ctx.strokeStyle = "blue";
    ctx.stroke();

    // Draw ball at current position
    let [ballX, ballZ] = traj[i];
    ctx.beginPath();
    ctx.arc(ballX * scaleX, canvas.height - ballZ * scaleY, 5, 0, Math.PI * 2);
    console.log("Ball position in pixels is: " + ballX * scaleX, canvas.height - ballZ * scaleY);
    console.log("Variables: " + ballX, scaleX, canvas.height, ballZ, scaleY);
    ctx.fillStyle = "red";
    ctx.fill();

    i++;
    if (i < traj.length) {
      requestAnimationFrame(drawFrame); // schedule next frame
    }
  }

  console.log("animating traj");

  drawFrame();
}

